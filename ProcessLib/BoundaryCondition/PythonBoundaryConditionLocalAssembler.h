/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "PythonBoundaryCondition.h"

#include "NumLib/DOF/DOFTableUtil.h"

#include "GenericNaturalBoundaryConditionLocalAssembler.h"
#include "PythonBoundaryConditionDetail.h"

namespace ProcessLib
{
template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class PythonBoundaryConditionLocalAssembler final
    : public GenericNaturalBoundaryConditionLocalAssembler<
          ShapeFunction, IntegrationMethod, GlobalDim>
{
    using Base = GenericNaturalBoundaryConditionLocalAssembler<
        ShapeFunction, IntegrationMethod, GlobalDim>;

public:
    PythonBoundaryConditionLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        bool is_axially_symmetric,
        unsigned const integration_order,
        PythonBoundaryConditionData const& data)
        : Base(e, is_axially_symmetric, integration_order),
          _data(data),
          _element(e)
    {
    }

    void assemble(std::size_t const boundary_element_id,
                  NumLib::LocalToGlobalIndexMap const& dof_table_boundary,
                  double const t, const GlobalVector& x, GlobalMatrix& /*K*/,
                  GlobalVector& b, GlobalMatrix* Jac) override
    {
        using ShapeMatricesType =
            ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
        using FemType =
            NumLib::TemplateIsoparametric<ShapeFunction, ShapeMatricesType>;
        FemType fe(*static_cast<const typename ShapeFunction::MeshElement*>(
            &_element));

        unsigned const num_integration_points =
            Base::_integration_method.getNumberOfPoints();
        auto const num_var = _data.dof_table_bulk.getNumberOfVariables();
        auto const num_nodes = _element.getNumberOfNodes();

        auto const& bulk_node_ids_map =
            *_data.boundary_mesh.getProperties().getPropertyVector<std::size_t>(
                "bulk_node_ids");

        // gather primary variables
        // TODO there might be problems with mixed ansatz functions
        std::vector<double> primary_variables;
        for (int var = 0; var < num_var; ++var)
        {
            auto const num_comp =
                _data.dof_table_bulk.getNumberOfVariableComponents(var);
            for (int comp = 0; comp < num_comp; ++comp)
            {
                for (unsigned element_node_id = 0; element_node_id < num_nodes;
                     ++element_node_id)
                {
                    auto* node = _element.getNode(element_node_id);
                    auto const boundary_node_id = node->getID();
                    auto const bulk_node_id =
                        bulk_node_ids_map[boundary_node_id];
                    MeshLib::Location loc{_data.bulk_mesh_id,
                                          MeshLib::MeshItemType::Node,
                                          bulk_node_id};
                    auto const dof_idx =
                        _data.dof_table_bulk.getGlobalIndex(loc, var, comp);
                    primary_variables.push_back(x[dof_idx]);
                }
            }
        }

        auto const num_comp_total =
            _data.dof_table_bulk.getNumberOfComponents();
        if (num_nodes * num_comp_total != primary_variables.size())
            OGS_FATAL("size mismatch");

        // TODO check if mapping is correct (row vs col major)
        Eigen::Map<Eigen::MatrixXd> primary_variables_mat(
            primary_variables.data(), num_nodes, num_comp_total);

        Eigen::VectorXd local_rhs = Eigen::VectorXd::Zero(num_nodes);
        Eigen::MatrixXd local_Jac =
            Eigen::MatrixXd::Zero(num_nodes, num_nodes * num_comp_total);
        bool has_dFlux = true;

        for (unsigned ip = 0; ip < num_integration_points; ip++)
        {
            auto const& sm = Base::_shape_matrices[ip];
            auto const coords = fe.interpolateCoordinates(sm.N);
            Eigen::VectorXd prim_vars =
                sm.N * primary_variables_mat;  // TODO problems with mixed
                                               // ansatz functions
            auto const res = _data.bc_object->getFlux(t, coords, prim_vars);
            if (!_data.bc_object->isOverriddenNatural())
                throw MethodNotOverriddenInDerivedClassException{};
            // TODO better documentation
            if (!std::get<0>(res))
                return;

            auto const& wp = Base::_integration_method.getWeightedPoint(ip);
            auto const w = sm.detJ * wp.getWeight() * sm.integralMeasure;
            local_rhs.noalias() += sm.N * (std::get<1>(res) * w);

            auto const& dFlux = std::get<2>(res);
            if (static_cast<int>(dFlux.size()) != num_comp_total)
            {
                // This strict check is technically mandatory only if a Jacobian
                // is assembled. However, it is done as a consistency check also
                // for cases without Jacobian assembly.
                OGS_FATAL(
                    "The Python BC must return the derivative of the flux "
                    "w.r.t. each primary variable. %d components expected. %d "
                    "components returned from Python.",
                    num_comp_total, dFlux.size());
            }
            if (Jac)
            {
                for (int comp = 0; comp < num_comp_total; ++comp)
                {
                    auto const top = 0;
                    auto const left = comp * num_nodes;
                    auto const width = num_nodes;
                    auto const height = num_nodes;
                    // The assignement -= takes into account the sign convention
                    // of 1st-order in time ODE systems in OpenGeoSys.
                    local_Jac.block(top, left, width, height).noalias() -=
                        sm.N.transpose() * (dFlux[comp] * w) * sm.N;
                }
            }
            else
            {
                has_dFlux = false;
            }
        }

        // TODO change for Newton (and Picard)!
        auto const& indices_comp =
            dof_table_boundary(boundary_element_id, _data.global_component_id)
                .rows;
        b.add(indices_comp, local_rhs);

        if (has_dFlux && Jac)
        {
            // only assemble a block of the Jacobian, not the whole local matrix
            auto const indices_all =
                NumLib::getIndices(boundary_element_id, dof_table_boundary);
            MathLib::RowColumnIndices<GlobalIndexType> rci{indices_comp,
                                                           indices_all};
            Jac->add(rci, local_Jac);
        }
    }

private:
    PythonBoundaryConditionData const& _data;
    MeshLib::Element const& _element;
};

}  // namespace ProcessLib
