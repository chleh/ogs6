/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_GENERICNATURALBOUNDARYCONDITIONLOCALASSEMBLER_H
#define PROCESSLIB_GENERICNATURALBOUNDARYCONDITIONLOCALASSEMBLER_H

#include "NumLib/Fem/ShapeMatrixPolicy.h"

#include "NumLib/DOF/DOFTableUtil.h"
#include "ProcessLib/Parameter/Parameter.h"

#include "ProcessLib/Utils/InitShapeMatrices.h"

namespace ProcessLib
{
class BalanceProcessLocalAssemblerInterface
{
public:
    virtual ~BalanceProcessLocalAssemblerInterface() = default;

    virtual void assemble(
        std::size_t const id,
        NumLib::LocalToGlobalIndexMap const& dof_table_boundary, double const t,
        const GlobalVector& x, GlobalMatrix& K, GlobalVector& b) = 0;
};

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class BalanceProcessLocalAssembler final
    : public BalanceProcessLocalAssemblerInterface
{
protected:
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using NodalMatrixType = typename ShapeMatricesType::NodalMatrixType;
    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;

public:
    /// The neumann_bc_value factor is directly integrated into the local
    /// element matrix.
    NeumannBoundaryConditionLocalAssembler(
        MeshLib::Element const& surface_element,
        std::size_t const local_matrix_size,
        unsigned const integration_order,
        ... // surface mesh, bulk mesh...
        )
        :
          // TODO: integration points for surface mesh and
          // shape matrices for bulk mesh
          _shape_matrices(initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                            IntegrationMethod, GlobalDim>(
              e /* bulk mesh element */,
                              integration_order)),
          _integration_order(integration_order),
          _neumann_bc_parameter(neumann_bc_parameter),
          _local_rhs(local_matrix_size)
    {
    }

    void assemble(std::size_t const id,
                  NumLib::LocalToGlobalIndexMap const& dof_table_boundary,
                  double const /*t*/, const GlobalVector& x,
                  GlobalMatrix& /*K*/, GlobalVector& /*b*/) override
    {
        // Ausgabe direkt in Property Vector vom Surface Mesh schreiben.

        // lokale Knotenwerte:
        auto const indices = NumLib::getIndices(mesh_item_id, dof_table);
        auto const local_x = x.get(indices);


        // Interpolation m.H. von Shapefunctions (mglw. nicht nötig)
        double p;
        NumLib::shapeFunctionInterpolate(local_x, _shape_matrices[0].N, p);

        //        for (std::size_t ip(0); ip < n_integration_points; ip++)
        //        {

        // Darcy velocity
        using GlobalDimVectorType =
            typename ShapeMatricesType::GlobalDimVectorType;
        GlobalDimVectorType const darcy_velocity =
            -k * sm.dNdx * Eigen::Map<const NodalVectorType>(
                               local_x.data(), ShapeFunction::NPOINTS);

        // } end for

//        IntegrationMethod integration_method(Base::_integration_order);
//        std::size_t const n_integration_points =
//            integration_method.getNumberOfPoints();

//        SpatialPosition pos;
//        pos.setElementID(id);

//        for (std::size_t ip(0); ip < n_integration_points; ip++)
//        {
//            pos.setIntegrationPoint(ip);
//            auto const& sm = Base::_shape_matrices[ip];
//            auto const& wp = integration_method.getWeightedPoint(ip);
//            _local_rhs.noalias() +=
//                sm.N * _neumann_bc_parameter.getTuple(t, pos).front() *
//                sm.detJ * wp.getWeight();
//        }

//        auto const indices = NumLib::getIndices(id, dof_table_boundary);
//        b.add(indices, _local_rhs);
    }

private:
    // vom bulk element an den surface-integrationspunkten
    std::vector<typename ShapeMatricesType::ShapeMatrices> const
        _shape_matrices;

    unsigned const _integration_order;
};

}  // ProcessLib

#endif  // PROCESSLIB_GENERICNATURALBOUNDARYCONDITIONLOCALASSEMBLER_H
