/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef PROCESS_LIB_TES_FEM_IMPL_H_
#define PROCESS_LIB_TES_FEM_IMPL_H_

#include "TESFEM.h"

#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"

// #include "logog/include/logog.hpp"

namespace ProcessLib
{

namespace TES
{


template <typename ShapeFunction_,
          typename IntegrationMethod_,
          typename GlobalMatrix,
          typename GlobalVector,
          unsigned GlobalDim>
void
LocalAssemblerData<ShapeFunction_,
    IntegrationMethod_,
    GlobalMatrix,
    GlobalVector,
    GlobalDim>::
init(MeshLib::Element const& e,
     std::size_t const /*local_matrix_size*/,
     unsigned const integration_order, const TESProcessInterface* process)
{
    using FemType = NumLib::TemplateIsoparametric<ShapeFunction, ShapeMatricesType>;

    FemType fe(*static_cast<const typename ShapeFunction::MeshElement*>(&e));


    _integration_order = integration_order;
    IntegrationMethod_ integration_method(_integration_order);
    std::size_t const n_integration_points = integration_method.getNPoints();

    _shape_matrices.resize(n_integration_points);
    for (std::size_t ip(0); ip < n_integration_points; ip++)
    {
        _shape_matrices[ip].resize(ShapeFunction::DIM, ShapeFunction::NPOINTS);
        fe.computeShapeFunctions(
                    integration_method.getWeightedPoint(ip).getCoords(),
                    _shape_matrices[ip]);
    }

    _localA.resize(MAT_SIZE, MAT_SIZE);
    _localRhs.resize(MAT_SIZE);
    // _localA.reset(new NodalMatrixType);
    // _localRhs.reset(new NodalVectorType);

    // DBUG("local matrix size: %i", local_matrix_size);

    _data._process = process;
}


template <typename ShapeFunction_,
          typename IntegrationMethod_,
          typename GlobalMatrix,
          typename GlobalVector,
          unsigned GlobalDim>
void
LocalAssemblerData<ShapeFunction_,
    IntegrationMethod_,
    GlobalMatrix,
    GlobalVector,
    GlobalDim>::
assemble(std::vector<double> const& localX,
         std::vector<double> const& localSecondaryVariables)
{

    _localA.setZero();
    _localRhs.setZero();

    IntegrationMethod_ integration_method(_integration_order);
    unsigned const n_integration_points = integration_method.getNPoints();

    _data.preEachAssemble(n_integration_points);

    for (std::size_t ip(0); ip < n_integration_points; ip++)
    {
        auto const& sm = _shape_matrices[ip];
        auto const& wp = integration_method.getWeightedPoint(ip);
        auto const weight = wp.getWeight();

        _data.assembleIntegrationPoint(&_localA, &_localRhs, localX, localSecondaryVariables,
                                       sm.N, sm.dNdx, sm.detJ, weight);
    }

    // first timestep:
    const Eigen::Map<const Eigen::VectorXd> oldX(localX.data(), localX.size());
    _data.postEachAssemble(&_localA, &_localRhs, oldX);
}


template <typename ShapeFunction_,
          typename IntegrationMethod_,
          typename GlobalMatrix,
          typename GlobalVector,
          unsigned GlobalDim>
void
LocalAssemblerData<ShapeFunction_,
    IntegrationMethod_,
    GlobalMatrix,
    GlobalVector,
    GlobalDim>::
addToGlobal(GlobalMatrix& A, GlobalVector& rhs,
            AssemblerLib::LocalToGlobalIndexMap::RowColumnIndices const& indices) const
{
    A.add(indices, _localA);
    rhs.add(indices.rows, _localRhs);
}


}   // namespace TES
}   // namespace ProcessLib



#endif // PROCESS_LIB_TES_FEM_IMPL_H_

