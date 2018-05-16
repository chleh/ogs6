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
class PythonConditionLocalAssembler final
    : public GenericNaturalBoundaryConditionLocalAssembler<
          ShapeFunction, IntegrationMethod, GlobalDim>
{
    using Base = GenericNaturalBoundaryConditionLocalAssembler<
        ShapeFunction, IntegrationMethod, GlobalDim>;

public:
    /// The neumann_bc_value factor is directly integrated into the local
    /// element matrix.
    PythonConditionLocalAssembler(MeshLib::Element const& e,
                                  std::size_t const local_matrix_size,
                                  bool is_axially_symmetric,
                                  unsigned const integration_order,
                                  PythonBoundaryConditionData const& data)
        : Base(e, is_axially_symmetric, integration_order),
          _data(data),
          _local_rhs(local_matrix_size),
          _element(e)
    {
    }

    void assemble(std::size_t const id,
                  NumLib::LocalToGlobalIndexMap const& dof_table_boundary,
                  double const t, const GlobalVector& x, GlobalMatrix& /*K*/,
                  GlobalVector& b) override
    {
        using ShapeMatricesType =
            ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
        using FemType =
            NumLib::TemplateIsoparametric<ShapeFunction, ShapeMatricesType>;
        FemType fe(*static_cast<const typename ShapeFunction::MeshElement*>(
            &_element));

        auto* bc =
            _data.scope[_data.bc_object.c_str()].cast<::PyBoundaryCondition*>();

        _local_rhs.setZero();

        unsigned const n_integration_points =
            Base::_integration_method.getNumberOfPoints();

        // std::vector<double> primary_variables;
        auto const num_var = _data.dof_table_boundary->getNumberOfVariables();
        auto const num_comp_total =
            _data.dof_table_boundary->getNumberOfComponents();
        DBUG("num var %i num comp total %i.", num_var, num_comp_total);

        auto const num_var2 = _data.dof_table_bulk.getNumberOfVariables();
        auto const num_comp_total2 =
            _data.dof_table_bulk.getNumberOfComponents();
        DBUG("num var2 %i num comp total2 %i.", num_var2, num_comp_total2);
        DBUG("num nodes %i.", _element.getNumberOfNodes());

        // Gathering primary variables
        // TODO there might be problems with mixed ansatz functions
        std::vector<double> primary_variables;
        auto const num_nodes = _element.getNumberOfNodes();
        for (int var = 0; var < num_var2; ++var)
        {
            auto const num_comp =
                _data.dof_table_bulk.getNumberOfVariableComponents(var);
            for (int comp = 0; comp < num_comp; ++comp)
            {
                for (unsigned element_node_id = 0; element_node_id < num_nodes;
                     ++element_node_id)
                {
                    auto* node = _element.getNode(element_node_id);
                    auto const node_id = node->getID();
                    MeshLib::Location loc{_data.mesh.getID(),
                                          MeshLib::MeshItemType::Node, node_id};
                    auto const dof_idx =
                        _data.dof_table_bulk.getGlobalIndex(loc, var, comp);
                    primary_variables.push_back(x[dof_idx]);
                }
            }
        }
        DBUG("prim var");
        for (auto v : primary_variables)
        {
            DBUG("  prim var %g.", v);
        }
        DBUG("prim var end");

        if (num_nodes * num_comp_total2 != primary_variables.size())
            OGS_FATAL("size mismatch");

        // TODO check if mapping is correct
        Eigen::Map<Eigen::MatrixXd> primary_variables_mat(
            primary_variables.data(), num_nodes, num_comp_total2);

        std::vector<GlobalIndexType> indices2 =
            NumLib::getIndices(_element.getID(), *_data.dof_table_boundary);

        DBUG("indices");
        for (auto i : indices2)
        {
            DBUG("index: %i.", i);
        }
        DBUG("end indices");

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            auto const& sm = Base::_shape_matrices[ip];
            auto const coords = fe.interpolateCoordinates(sm.N);
            Eigen::VectorXd prim_vars =
                sm.N * primary_variables_mat;  // TODO problems with mixed
                                               // ansatz functions
            auto const res = bc->getFlux(t, coords, prim_vars);
            if (!res.first)
                return;
            if (!bc->isOverriddenNatural())
                throw PyNotOverridden{};

            auto const& wp = Base::_integration_method.getWeightedPoint(ip);
            _local_rhs.noalias() += sm.N * res.second * sm.detJ *
                                    wp.getWeight() * sm.integralMeasure;
        }

        auto const indices = NumLib::getIndices(id, dof_table_boundary);
        DBUG("id %i vs. %i.", id, _element.getID());
        DBUG("indices size %i vs. %i.", indices.size(), indices2.size());
        b.add(indices, _local_rhs);
    }

private:
    PythonBoundaryConditionData const& _data;
    typename Base::NodalVectorType _local_rhs;
    MeshLib::Element const& _element;

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace ProcessLib
