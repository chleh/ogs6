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
#include "NumLib/Function/Interpolation.h"

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
     unsigned const integration_order, TESProcessInterface* process)
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

    _data._AP = & process->getAssemblyParams();

    _data.init(n_integration_points, GlobalDim);
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
         std::vector<double> const& localXPrevTs)
{

    _localA.setZero();
    _localRhs.setZero();

    IntegrationMethod_ integration_method(_integration_order);
    unsigned const n_integration_points = integration_method.getNPoints();

    _data.preEachAssemble();

    for (std::size_t ip(0); ip < n_integration_points; ip++)
    {
        auto const& sm = _shape_matrices[ip];
        auto const& wp = integration_method.getWeightedPoint(ip);
        auto const weight = wp.getWeight();

        _data.assembleIntegrationPoint(ip, &_localA, &_localRhs, localX,
                                       sm.N, sm.dNdx, sm.J, sm.detJ, weight);
    }

    // first timestep:
    const Eigen::Map<const Eigen::VectorXd> oldX(localXPrevTs.data(), localXPrevTs.size());
    _data.postEachAssemble(_localA, _localRhs, oldX);
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


template <typename ShapeFunction_,
          typename IntegrationMethod_,
          typename GlobalMatrix,
          typename GlobalVector,
          unsigned GlobalDim>
std::shared_ptr<const std::vector<double> >
LocalAssemblerData<ShapeFunction_,
    IntegrationMethod_,
    GlobalMatrix,
    GlobalVector,
    GlobalDim>::
getIntegrationPointValues(SecondaryVariables var, NumLib::LocalNodalDOF& nodal_dof) const
{

    switch (var)
    {
    case SecondaryVariables::REACTION_RATE:
    case SecondaryVariables::SOLID_DENSITY:
    case SecondaryVariables::VELOCITY_X:
    case SecondaryVariables::VELOCITY_Y:
    case SecondaryVariables::VELOCITY_Z:
    case SecondaryVariables::REACTION_KINETIC_INDICATOR:
    case SecondaryVariables::LOADING:
        // These cases do not need access to nodal values
        // Thus, they can be handled inside _data
        return _data.getIntegrationPointValues(var);

    // TODO [CL] the following cases could be better provided directly using nodal values without extrapolation
    case SecondaryVariables::VAPOUR_PARTIAL_PRESSURE:
    {
        IntegrationMethod_ integration_method(_integration_order);
        auto const n_integration_points = integration_method.getNPoints();

        auto pVs = std::make_shared<std::vector<double> >();
        pVs->reserve(n_integration_points);

        auto const ps = nodal_dof.getElementNodalValues(0); // TODO [CL] use constants for DOF indices
        auto const xs = nodal_dof.getElementNodalValues(2);

        auto const& AP = *_data._AP;

        for (auto const& sm : _shape_matrices)
        {
            double p, xm;

            using Array = std::array<double*, 1>;
            NumLib::shapeFunctionInterpolate(ps, sm.N, Array{ &p  });
            NumLib::shapeFunctionInterpolate(xs, sm.N, Array{ &xm });

            xm = Trafo::x(xm);

            auto const xn = AP._adsorption->get_molar_fraction(xm, AP._M_react, AP._M_inert);
            pVs->push_back(p * xn);
        }

        return pVs;
    }
    case SecondaryVariables::RELATIVE_HUMIDITY:
    {
        IntegrationMethod_ integration_method(_integration_order);
        auto const n_integration_points = integration_method.getNPoints();

        auto rhs = std::make_shared<std::vector<double> >();
        rhs->reserve(n_integration_points);

        auto const nodal_vals = nodal_dof.getElementNodalValues();

        auto const& AP = *_data._AP;

        for (auto const& sm : _shape_matrices)
        {
            double p, T, xm;

            using Array = std::array<double*, 3>;
            NumLib::shapeFunctionInterpolate(nodal_vals, sm.N, Array{ &p, &T, &xm });

            xm = Trafo::x(xm);

            auto const xn = AP._adsorption->get_molar_fraction(xm, AP._M_react, AP._M_inert);
            auto const pS = AP._adsorption->get_equilibrium_vapour_pressure(T);
            rhs->push_back(p * xn / pS);
        }

        return rhs;
    }
    case SecondaryVariables::EQUILIBRIUM_LOADING:
    {
        IntegrationMethod_ integration_method(_integration_order);
        auto const n_integration_points = integration_method.getNPoints();

        auto Cs = std::make_shared<std::vector<double> >();
        Cs->reserve(n_integration_points);

        auto const nodal_vals = nodal_dof.getElementNodalValues();

        auto const& AP = *_data._AP;

        for (auto const& sm : _shape_matrices)
        {
            double p, T, xm;

            using Array = std::array<double*, 3>;
            NumLib::shapeFunctionInterpolate(nodal_vals, sm.N, Array{ &p, &T, &xm });

            xm = Trafo::x(xm);

            auto const xn = AP._adsorption->get_molar_fraction(xm, AP._M_react, AP._M_inert);
            auto const pV = p * xn;
            if (pV < 0.0) {
                Cs->push_back(0.0);
            } else {
                Cs->push_back(AP._adsorption->get_equilibrium_loading(pV, T, AP._M_react));
            }
        }

        return Cs;
    }
    case SecondaryVariables::REACTION_DAMPING_FACTOR:
    {
        auto alphas = std::make_shared<std::vector<double> >();
        alphas->resize(_shape_matrices.size(), _data.reaction_damping_factor);

        return alphas;
    }
    }


    return nullptr;
}


template <typename ShapeFunction_,
          typename IntegrationMethod_,
          typename GlobalMatrix,
          typename GlobalVector,
          unsigned GlobalDim>
bool
LocalAssemblerData<ShapeFunction_,
    IntegrationMethod_,
    GlobalMatrix,
    GlobalVector,
    GlobalDim>::
checkBounds(std::vector<double> const& localX,
            const std::vector<double>& localX_pts)
{
    double alpha = 1.0;

    const double min_xmV = 1e-6;
    const std::size_t nnodes = localX.size() / NODAL_DOF;
    const std::size_t xmV_offset = (NODAL_DOF - 1)*nnodes;

    for (std::size_t i=0; i<nnodes; ++i)
    {
        auto const xnew = localX[xmV_offset+i];
        if (xnew < min_xmV)
        {
            auto const xold = localX_pts[i+xmV_offset];
            const auto a = xold / (xold - xnew);
            if (a<alpha) DBUG("xo %g, xn %g, a %g", xold, xnew, a);
            alpha = std::min(alpha, a);
            _data.bounds_violation[i] = true;
        }
        else if (xnew > 1.0)
        {
            auto const xold = localX_pts[i+xmV_offset];
            const auto a = xold / (xnew - xold);
            if (a<alpha) DBUG("xo %g, xn %g, a %g", xold, xnew, a);
            alpha = std::min(alpha, a);
            _data.bounds_violation[i] = true;
        }
        else
        {
            _data.bounds_violation[i] = false;
        }
    }

    assert (alpha > 0.0);

    if (alpha != 1.0)
    {
        if (_data._AP->_number_of_try_of_iteration <=2) {
            _data.reaction_damping_factor *= sqrt(std::min(alpha, 0.5));
                                            // * sqrt(_data.reaction_damping_factor);
        } else {
            _data.reaction_damping_factor *= std::min(alpha, 0.5);
        }
    }

    DBUG("new damping factor: %g", _data.reaction_damping_factor);

    return alpha == 1.0;
}


}   // namespace TES
}   // namespace ProcessLib



#endif // PROCESS_LIB_TES_FEM_IMPL_H_

