/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * The code of this file is used to decouple the evaluation of matrix elements
 * from the rest of OGS6,
 * not all of OGS6 has to be recompiled every time a small change is done.
 */

#pragma once

#include <logog/include/logog.hpp>

#include "MaterialLib/Adsorption/Adsorption.h"
#include "MathLib/Nonlinear/Root1D.h"
#include "NumLib/Function/Interpolation.h"

#include "TESDielectricByKraus.h"
#include "TESLocalAssemblerInner-fwd.h"
#include "TESOGS5MaterialModels.h"

namespace ProcessLib
{
namespace TES
{
template <typename Traits>
TESLocalAssemblerInner<Traits>::TESLocalAssemblerInner(
    const AssemblyParams& ap, const std::size_t element_id,
    const unsigned num_int_pts, const unsigned dimension)
    : _d{ap, element_id, num_int_pts, dimension}
{
}

template <typename Traits>
Eigen::Matrix3d TESLocalAssemblerInner<Traits>::getMassCoeffMatrix(
    const unsigned int_pt)
{
    // TODO: Dalton's law property
    const double dxn_dxm = Adsorption::AdsorptionReaction::dMolarFraction(
        _d.vapour_mass_fraction, _d.ap.M_react, _d.ap.M_inert);
    const double R = MaterialLib::PhysicalConstant::IdealGasConstant;

    const double M_pp = _d.poro / _d.p * _d.rho_GR;
    const double M_pT = -_d.poro / _d.T * _d.rho_GR;
    const double M_px =
        (_d.ap.M_react - _d.ap.M_inert) * _d.p / (R * _d.T) * dxn_dxm * _d.poro;

    const double M_Tp = -_d.poro;
    const double M_TT =
        _d.poro * _d.rho_GR * _d.ap.cpG  // TODO: vapour heat capacity
        + (1.0 - _d.poro) *
              _d.ap.reactive_solid->getSolidDensity(
                  *_d.reactive_solid_state[int_pt]) *
              _d.ap.cpS;  // TODO: adsorbate heat capacity
    const double M_Tx = 0.0;

    const double M_xp = 0.0;
    const double M_xT = 0.0;
    const double M_xx = _d.poro * _d.rho_GR;

    Eigen::Matrix3d M;
    M << M_pp, M_pT, M_px, M_Tp, M_TT, M_Tx, M_xp, M_xT, M_xx;

    return M;
}

template <typename Traits>
typename Traits::LaplaceMatrix
TESLocalAssemblerInner<Traits>::getLaplaceCoeffMatrix(const double t,
                                                      const unsigned int_pt,
                                                      const unsigned dim)
{
    const double eta_GR = fluid_viscosity(_d.p, _d.T, _d.vapour_mass_fraction);

    const double lambda_F =
        fluid_heat_conductivity(_d.p, _d.T, _d.vapour_mass_fraction);
    const double lambda_S = _d.ap.solid_heat_cond;

    using Mat = typename Traits::MatrixDimDim;

    typename Traits::LaplaceMatrix L =
        Traits::LaplaceMatrix::Zero(dim * NODAL_DOF, dim * NODAL_DOF);

    // TODO: k_rel
    // L_pp
    SpatialPosition pos;
    pos.setElementID(_d.element_id);
    pos.setIntegrationPoint(int_pt);
    auto const& perm = *_d.ap.permeability;
    auto const k = perm(t, pos).front();
    Traits::blockDimDim(L, 0, 0, dim, dim) =
        Mat::Identity(dim, dim) * (k * _d.rho_GR / eta_GR);

    // TODO: add zeolite part
    // L_TT
    Traits::blockDimDim(L, dim, dim, dim, dim) =
        Mat::Identity(dim, dim) *
        (_d.poro * lambda_F + (1.0 - _d.poro) * lambda_S);

    // L_xx
    auto const D =
        _d.ap.diffusion_coefficient_component->getDiffusionCoefficient(_d.p,
                                                                       _d.T);

    Traits::blockDimDim(L, 2 * dim, 2 * dim, dim, dim) =
        Mat::Identity(dim, dim) * (_d.ap.tortuosity * _d.poro * _d.rho_GR * D);

    return L;
}

template <typename Traits>
Eigen::Matrix3d TESLocalAssemblerInner<Traits>::getAdvectionCoeffMatrix(
    const unsigned /*int_pt*/)
{
    const double A_pp = 0.0;
    const double A_pT = 0.0;

    const double A_px = 0.0;

    const double A_Tp = 0.0;

    const double A_TT = _d.rho_GR * _d.ap.cpG;  // porosity?
    const double A_Tx = 0.0;

    const double A_xp = 0.0;
    const double A_xT = 0.0;
    const double A_xx = _d.rho_GR;  // porosity?

    Eigen::Matrix3d A;
    A << A_pp, A_pT, A_px, A_Tp, A_TT, A_Tx, A_xp, A_xT, A_xx;

    return A;
}

template <typename Traits>
Eigen::Matrix3d TESLocalAssemblerInner<Traits>::getContentCoeffMatrix(
    const unsigned int_pt)
{
    const double C_pp = 0.0;
    const double C_pT = 0.0;

    const double C_px = 0.0;

    const double C_Tp = 0.0;

    const double C_TT = 0.0;
    const double C_Tx = 0.0;

    const double C_xp = 0.0;
    const double C_xT = 0.0;
    const double C_xx = (_d.poro - 1.0) * _d.ap.reactive_solid->getOverallRate(
                                              *_d.reaction_rate[int_pt]);

    Eigen::Matrix3d C;
    C << C_pp, C_pT, C_px, C_Tp, C_TT, C_Tx, C_xp, C_xT, C_xx;

    return C;
}

template <typename Traits>
Eigen::Vector3d TESLocalAssemblerInner<Traits>::getRHSCoeffVector(
    const unsigned int_pt)
{
    auto const reaction_enthalpy = _d.ap.reaction_rate->getHeatOfReaction(
        _d.p_V, _d.T, _d.reactive_solid_state[int_pt].get());

    const double rhs_p =
        (_d.poro - 1.0) *
        _d.ap.reactive_solid->getOverallRate(
            *_d.reaction_rate[int_pt]);  // TODO [CL] body force term

    double rhs_T =
        _d.rho_GR * _d.poro * _d.ap.fluid_specific_heat_source +
        (1.0 - _d.poro) * _d.ap.reactive_solid->getHeatingRate(
                              reaction_enthalpy, *_d.reaction_rate[int_pt]) +
        _d.ap.reactive_solid->getSolidDensity(
            *_d.reactive_solid_state[int_pt]) *
            (1.0 - _d.poro) * _d.ap.solid_specific_heat_source;

    if (_d.ap.volumetric_heat_loss)
        rhs_T += _d.ap.volumetric_heat_loss->getHeatLoss(_d.T);

    // TODO [CL] momentum production term

    if (_d.ap.dielectric_heating_term_enabled)
    {
        auto const loading = Adsorption::AdsorptionReaction::getLoading(
            _d.ap.reactive_solid->getSolidDensity(
                *_d.reactive_solid_state[int_pt]),
            _d.ap.rho_SR_dry);

        rhs_T += _d.ap.heating_power_scaling.getValue(_d.ap.current_time) *
                 (1.0 - _d.poro) *
                 getVolumetricJouleHeatingPower(_d.T, loading);
    }

    const double rhs_x = (_d.poro - 1.0) * _d.ap.reactive_solid->getOverallRate(
                                               *_d.reaction_rate[int_pt]);

    Eigen::Vector3d rhs;
    rhs << rhs_p, rhs_T, rhs_x;

    return rhs;
}

template <typename Traits>
void TESLocalAssemblerInner<Traits>::preEachAssembleIntegrationPoint(
    const unsigned int_pt,
    const double t,
    const std::vector<double>& localX,
    typename Traits::ShapeMatrices const& sm)
{
#ifndef NDEBUG
    // fill local data with garbage to aid in debugging
    _d.p = _d.T = _d.vapour_mass_fraction = _d.p_V = _d.rho_GR = _d.poro =
        std::numeric_limits<double>::quiet_NaN();
#endif

    SpatialPosition pos;
    pos.setElementID(_d.element_id);
    pos.setIntegrationPoint(int_pt);
    _d.poro = (*_d.ap.poro)(t, pos).front();

    NumLib::shapeFunctionInterpolate(localX, sm.N, _d.p, _d.T,
                                     _d.vapour_mass_fraction);

    _d.p_V = _d.p * Adsorption::AdsorptionReaction::getMolarFraction(
                        _d.vapour_mass_fraction, _d.ap.M_react, _d.ap.M_inert);

    if (_d.ap.reaction_rate->computeReactionRate(
            _d.ap.delta_t, _d.p, _d.T, _d.p_V, *_d.reaction_rate[int_pt],
            _d.reactive_solid_state[int_pt].get(),
            _d.reaction_rate_data[int_pt].get()))
    {
        _d.vapour_mass_fraction =
            Adsorption::AdsorptionReaction::getMassFraction(
                _d.p_V / _d.p, _d.ap.M_react, _d.ap.M_inert);
    }

    assert(_d.p > 0.0);
    assert(_d.T > 0.0);
    assert(0.0 <= _d.vapour_mass_fraction && _d.vapour_mass_fraction <= 1.0);

    _d.rho_GR = fluid_density(_d.p, _d.T, _d.vapour_mass_fraction);
}

template <typename Traits>
void TESLocalAssemblerInner<Traits>::assembleIntegrationPoint(
    unsigned integration_point,
    const double t,
    std::vector<double> const& localX,
    typename Traits::ShapeMatrices const& sm,
    const double weight,
    Eigen::Map<typename Traits::LocalMatrix>& local_M,
    Eigen::Map<typename Traits::LocalMatrix>& local_K,
    Eigen::Map<typename Traits::LocalVector>& local_b)
{
    preEachAssembleIntegrationPoint(integration_point, t, localX, sm);

    auto const N =
        static_cast<unsigned>(sm.dNdx.cols());  // number of integration points
    auto const D =
        static_cast<unsigned>(sm.dNdx.rows());  // global dimension: 1, 2 or 3

    assert(N * NODAL_DOF == local_M.cols());

    auto const laplaceCoeffMat = getLaplaceCoeffMatrix(t, integration_point, D);
    assert(laplaceCoeffMat.cols() == D * NODAL_DOF);
    auto const massCoeffMat = getMassCoeffMatrix(integration_point);
    auto const advCoeffMat = getAdvectionCoeffMatrix(integration_point);
    auto const contentCoeffMat = getContentCoeffMatrix(integration_point);

    auto const velocity =
        (Traits::blockDimDim(laplaceCoeffMat, 0, 0, D, D) *
         (sm.dNdx *
          Eigen::Map<const typename Traits::Vector1Comp>(localX.data(),
                                                         N)  // grad_p
          / -_d.rho_GR))
            .eval();
    assert(velocity.size() == D);

    auto const detJ_w_im_NT =
        (sm.detJ * weight * sm.integralMeasure * sm.N.transpose()).eval();
    auto const detJ_w_im_NT_N = (detJ_w_im_NT * sm.N).eval();
    assert(detJ_w_im_NT_N.rows() == N && detJ_w_im_NT_N.cols() == N);

    auto const detJ_w_im_NT_vT_dNdx =
        (detJ_w_im_NT * velocity.transpose() * sm.dNdx).eval();
    assert(detJ_w_im_NT_vT_dNdx.rows() == N &&
           detJ_w_im_NT_vT_dNdx.cols() == N);

    for (unsigned r = 0; r < NODAL_DOF; ++r)
    {
        for (unsigned c = 0; c < NODAL_DOF; ++c)
        {
            Traits::blockShpShp(local_K, N * r, N * c, N, N).noalias() +=
                sm.detJ * weight * sm.integralMeasure * sm.dNdx.transpose() *
                    Traits::blockDimDim(laplaceCoeffMat, D * r, D * c, D, D) *
                    sm.dNdx  // end Laplacian part
                + detJ_w_im_NT_N * contentCoeffMat(r, c) +
                detJ_w_im_NT_vT_dNdx * advCoeffMat(r, c);
            Traits::blockShpShp(local_M, N * r, N * c, N, N).noalias() +=
                detJ_w_im_NT_N * massCoeffMat(r, c);
        }
    }

    auto const rhsCoeffVector = getRHSCoeffVector(integration_point);

    for (unsigned r = 0; r < NODAL_DOF; ++r)
    {
        Traits::blockShp(local_b, N * r, N).noalias() +=
            rhsCoeffVector(r) * sm.N.transpose() * sm.detJ * weight *
            sm.integralMeasure;
    }
}

template <typename Traits>
void TESLocalAssemblerInner<Traits>::preTimestep()
{
    // TODO check
    /*
    std::vector<double> solid_density;
    solid_density.reserve(_d.reactive_solid_state.size());
    for (auto& s : _d.reactive_solid_state) {
        solid_density.emplace_back(s->solid_density);
    }
    */

    for (std::size_t i = 0; i < _d.reactive_solid_state.size(); ++i)
    {
        if (_d.reactive_solid_state[i])
            _d.ap.reactive_solid->preTimestep(*_d.reactive_solid_state[i],
                                              *_d.reaction_rate[i]);

        if (_d.reaction_rate_data[i])
            _d.reaction_rate_data[i]->preTimestep(*_d.reactive_solid_state[i],
                                                  *_d.reaction_rate[i]);
    }
}

}  // namespace TES
}  // namespace ProcessLib
