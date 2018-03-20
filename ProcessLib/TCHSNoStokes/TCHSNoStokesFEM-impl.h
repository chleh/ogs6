/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   TCHSNoStokesFEM-impl.h
 *  Created on November 29, 2017, 2:03 PM
 */

#pragma once

#include "TCHSNoStokesFEM.h"

#include <boost/math/special_functions/pow.hpp>

#include "MaterialLib/Adsorption/Adsorption.h"
#include "MathLib/KelvinVector.h"
#include "NumLib/Fem/CoordinatesMapping/NaturalNodeCoordinates.h"
#include "NumLib/Function/Interpolation.h"
#include "ProcessLib/CoupledSolutionsForStaggeredScheme.h"

#include "ProcessLib/TES/TESOGS5MaterialModels.h"

namespace ProcessLib
{
namespace TCHSNoStokes
{
template <typename ShapeFunctionVelocity, typename ShapeFunctionPressure,
          typename IntegrationMethod, int VelocityDim>
TCHSNoStokesLocalAssembler<ShapeFunctionVelocity, ShapeFunctionPressure,
                           IntegrationMethod, VelocityDim>::
    TCHSNoStokesLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        TCHSNoStokesProcessData<VelocityDim>& process_data)
    : _process_data(process_data),
      _integration_method(integration_order),
      _element(e),
      _is_axially_symmetric(is_axially_symmetric)
{
    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    _ip_data.resize(n_integration_points);

    auto const shape_matrices_linear =
        initShapeMatrices<ShapeFunctionPressure, ShapeMatricesTypePressure,
                          IntegrationMethod, VelocityDim>(
            e, is_axially_symmetric, _integration_method);

    SpatialPosition pos;
    pos.setElementID(_element.getID());

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        // velocity (subscript u)
        auto& ip_data = _ip_data[ip];
        auto const& sm_1 = shape_matrices_linear[ip];
        ip_data.integration_weight =
            _integration_method.getWeightedPoint(ip).getWeight() *
            sm_1.integralMeasure * sm_1.detJ;

        ip_data.N_1 = shape_matrices_linear[ip].N;
        ip_data.dNdx_1 = shape_matrices_linear[ip].dNdx;

        pos.setIntegrationPoint(ip);

        auto const mat_id = _process_data.material_ids[_element.getID()];
        auto& mat = _process_data.materials.at(mat_id);

        // TODO warning: 0.0 is the time!
        ip_data.reactive_solid_state =
            mat.reactive_solid->createReactiveSolidState(0.0 /* time */, pos);

        if (!mat.reaction_rate->isStateCompatible(
                *ip_data.reactive_solid_state))
            OGS_FATAL(
                "reaction rate and reactive solid state are incompatible.");

        ip_data.reaction_rate = mat.reactive_solid->createReactiveSolidRate();
        ip_data.reaction_rate_data =
            mat.reaction_rate->createReactionRateData();
    }
}

template <typename ShapeFunctionVelocity, typename ShapeFunctionPressure,
          typename IntegrationMethod, int VelocityDim>
void TCHSNoStokesLocalAssembler<
    ShapeFunctionVelocity, ShapeFunctionPressure, IntegrationMethod,
    VelocityDim>::assemble(double const t, std::vector<double> const& local_x,
                           std::vector<double>& local_M_data,
                           std::vector<double>& local_K_data,
                           std::vector<double>& local_rhs_data)
{
    assert(local_x.size() == Block::index(Block::X) + Block::size(Block::X));

    // make Eigen::Map<>s
    auto const nodal_p = Block::mapVectorSegment(local_x, Block::P);
    auto const nodal_T = Block::mapVectorSegment(local_x, Block::T);
    auto const nodal_xmV = Block::mapVectorSegment(local_x, Block::X);

    auto local_M = MathLib::createZeroedMatrix<
        typename ShapeMatricesTypeVelocity::template MatrixType<
            Block::index(Block::X) + Block::size(Block::X),
            Block::index(Block::X) + Block::size(Block::X)>>(
        local_M_data, Block::index(Block::X) + Block::size(Block::X),
        Block::index(Block::X) + Block::size(Block::X));

    auto local_K = MathLib::createZeroedMatrix<
        typename ShapeMatricesTypeVelocity::template MatrixType<
            Block::index(Block::X) + Block::size(Block::X),
            Block::index(Block::X) + Block::size(Block::X)>>(
        local_K_data, Block::index(Block::X) + Block::size(Block::X),
        Block::index(Block::X) + Block::size(Block::X));

    auto local_rhs = MathLib::createZeroedVector<
        typename ShapeMatricesTypeVelocity::template VectorType<
            Block::index(Block::X) + Block::size(Block::X)>>(
        local_rhs_data, Block::index(Block::X) + Block::size(Block::X));

    SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto const& ip_data = _ip_data[ip];

        // shape functions and weights /////////////////////////////////////////
        x_position.setIntegrationPoint(ip);
        auto const& w = ip_data.integration_weight;

        auto const& N_1 = ip_data.N_1;
        auto const& dNdx_1 = ip_data.dNdx_1;

        auto const x_coord =
            interpolateXCoordinate<ShapeFunctionVelocity,
                                   ShapeMatricesTypeVelocity>(_element, N_1);

        // interpolate nodal values ////////////////////////////////////////////
        // TODO fix
        Eigen::Matrix<double, VelocityDim, 1> const v_Darcy =
            Eigen::Matrix<double, VelocityDim, 1>::Zero();
        double const p = N_1 * nodal_p;
        double const T = N_1 * nodal_T;
        double x_mV = N_1 * nodal_xmV;

        // some shortcuts //////////////////////////////////////////////////////
        auto const N_1_T_N_1 = (N_1.transpose() * N_1).eval();
        auto const N_1_T_v_T_dNdx_1 =
            (N_1.transpose() * v_Darcy.transpose() * dNdx_1).eval();

        // material parameters /////////////////////////////////////////////////
        auto const mat_id = _process_data.material_ids[_element.getID()];
        auto const& mat = _process_data.materials.at(mat_id);

        double const M_R = mat.molar_mass_reactive;
        double const M_I = mat.molar_mass_inert;
        double p_V = p * Adsorption::AdsorptionReaction::getMolarFraction(
                             x_mV, M_R, M_I);

        // porosity
        auto const porosity = mat.porosity->getPorosity(x_coord);
        Eigen::Matrix<double, VelocityDim, 1> grad_porosity =
            Eigen::Matrix<double, VelocityDim, 1>::Zero();
        grad_porosity[0] = mat.porosity->getDPorosityDr(x_coord);

        // reaction
        if (mat.reaction_rate->computeReactionRate(
                _process_data.delta_t, p, T, p_V, *ip_data.reaction_rate,
                ip_data.reactive_solid_state.get(),
                ip_data.reaction_rate_data.get()))
        {
            x_mV = Adsorption::AdsorptionReaction::getMassFraction(p_V / p, M_R,
                                                                   M_I);
        }

        double const hat_rho_S =
            (1.0 - porosity) *
            mat.reactive_solid->getOverallRate(*ip_data.reaction_rate);
        auto const reaction_enthalpy = mat.reaction_rate->getHeatOfReaction(
            p_V, T, ip_data.reactive_solid_state.get());
        double const heating_rate =
            (1.0 - porosity) * mat.reactive_solid->getHeatingRate(
                                   reaction_enthalpy, *ip_data.reaction_rate);

        // fluid
        double const rho_GR = mat.fluid_density->getDensity(p, T, x_mV);
        double const alpha_T =
            mat.fluid_density->getThermalExpansionCoefficient(p, T, x_mV);
        double const beta_p = mat.fluid_density->getCompressibility(p, T, x_mV);
        double const gamma_x =
            mat.fluid_density->getDensityChangeWithComposition(p, T, x_mV);

        auto const mass_dispersion = mat.mass_dispersion->getMassDispersion(
            t, p, T, v_Darcy.norm(), x_coord, porosity,
            _process_data.probed_pressure, _process_data.probed_temperature,
            _process_data.probed_velocity);
        double const c_pG = mat.fluid_heat_capacity->getHeatCapacity(T, x_mV);

        // fluid viscosity/friction
        double const mu = mat.fluid_viscosity->getViscosity(p, T, x_mV);
        double const Re_0 =
            mat.reynolds_number->getRe(t, rho_GR, v_Darcy.norm(), mu);

        // solid
        double const rho_SR =
            mat.reactive_solid->getSolidDensity(*ip_data.reactive_solid_state);

        // fluid and solid
        double const total_heat_capacity =
            porosity * c_pG * rho_GR +
            (1.0 - porosity) *
                mat.solid_heat_capacity->getSpecificHeatCapacity(rho_SR, T) *
                rho_SR;
        auto const total_heat_conductivity =
            mat.heat_conductivity->getHeatConductivity(
                t, p, T, x_mV, x_coord, porosity, rho_GR, c_pG, Re_0,
                v_Darcy.norm(), _process_data.probed_velocity);

        // assemble local matrices /////////////////////////////////////////////

        // M_p? (total mass balance) ///////////////////////////////////////////
        // M_pp
        Block::block(local_M, Block::P, Block::P).noalias() +=
            N_1_T_N_1 * (rho_GR * porosity * beta_p * w);

        // M_pT
        Block::block(local_M, Block::P, Block::T).noalias() -=
            N_1_T_N_1 * (rho_GR * porosity * alpha_T * w);

        // M_px
        Block::block(local_M, Block::P, Block::X).noalias() +=
            N_1_T_N_1 * (rho_GR * porosity * gamma_x * w);

        // M_T? (total energy balance) /////////////////////////////////////////
        // M_Tp
        Block::block(local_M, Block::T, Block::P).noalias() -=
            N_1_T_N_1 * (porosity * T * alpha_T * w);

        // M_TT
        Block::block(local_M, Block::T, Block::T).noalias() +=
            N_1_T_N_1 * (total_heat_capacity * w);

        // M_Tx = 0

        // M_x? (vapour mass balance) //////////////////////////////////////////
        // M_xp = 0

        // M_xT = 0

        // M_xx
        Block::block(local_M, Block::X, Block::X).noalias() +=
            N_1_T_N_1 * (rho_GR * porosity * w);

        // M_v? (gas momentum balance) /////////////////////////////////////////
        // M_vp = 0

        // M_vT = 0

        // M_vx = 0

        // M_vv = 0

        // K_p? (total mass balance) ///////////////////////////////////////////
        // TODO check K_p?
        // K_pp
        Block::block(local_K, Block::P, Block::P).noalias() +=
            N_1_T_v_T_dNdx_1 * (rho_GR * beta_p * w);

        // K_pT
        Block::block(local_K, Block::P, Block::T).noalias() -=
            N_1_T_v_T_dNdx_1 * (rho_GR * alpha_T * w);

        // K_px
        Block::block(local_K, Block::P, Block::X).noalias() +=
            N_1_T_v_T_dNdx_1 * (rho_GR * gamma_x * w);

        // K_T? (total energy balance) /////////////////////////////////////////
        // K_Tp
        Block::block(local_K, Block::T, Block::T).noalias() +=
            N_1_T_v_T_dNdx_1 * ((1.0 - T * alpha_T) * w);

        // K_TT
        Block::block(local_K, Block::T, Block::T).noalias() +=
            N_1_T_v_T_dNdx_1 * (rho_GR * c_pG * w) +
            dNdx_1.transpose() * total_heat_conductivity * w * dNdx_1;

        // K_Tx = 0

        // K_x? (vapour mass balance) //////////////////////////////////////////
        // K_xp = 0

        // K_xT = 0

        // K_xx
        Block::block(local_K, Block::X, Block::X).noalias() +=
            N_1_T_v_T_dNdx_1 * (rho_GR * w) - N_1_T_N_1 * (hat_rho_S * w) +
            dNdx_1.transpose() * (rho_GR * w) * mass_dispersion * dNdx_1;

        // rhs /////////////////////////////////////////////////////////////////
        // rhs_p
        Block::segment(local_rhs, Block::P).noalias() -=
            N_1.transpose() * (hat_rho_S * w);

        // rhs_T
        Block::segment(local_rhs, Block::T).noalias() +=
            N_1.transpose() * (heating_rate * w);

        // rhs_x
        Block::segment(local_rhs, Block::X).noalias() -=
            N_1.transpose() * (hat_rho_S * w);
    }
}

template <typename ShapeFunctionVelocity, typename ShapeFunctionPressure,
          typename IntegrationMethod, int VelocityDim>
void TCHSNoStokesLocalAssembler<ShapeFunctionVelocity, ShapeFunctionPressure,
                                IntegrationMethod, VelocityDim>::
    postNonLinearSolverConcrete(std::vector<double> const&
                                /*local_x*/,
                                double const /*t*/,
                                bool const /*use_monolithic_scheme*/)
{
}

template <typename ShapeFunctionVelocity, typename ShapeFunctionPressure,
          typename IntegrationMethod, int VelocityDim>
std::vector<double> const&
TCHSNoStokesLocalAssembler<ShapeFunctionVelocity, ShapeFunctionPressure,
                           IntegrationMethod, VelocityDim>::
    getIntPtSolidDensity(const double /*t*/,
                         GlobalVector const& /*current_solution*/,
                         NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
                         std::vector<double>& cache) const
{
    cache.clear();
    cache.reserve(_ip_data.size());

    auto const mat_id = _process_data.material_ids[_element.getID()];
    auto const& mat = _process_data.materials.at(mat_id);

    for (auto& ip_data : _ip_data)
    {
        cache.emplace_back(
            mat.reactive_solid->getSolidDensity(*ip_data.reactive_solid_state));
    }
    return cache;
}

template <typename ShapeFunctionVelocity, typename ShapeFunctionPressure,
          typename IntegrationMethod, int VelocityDim>
std::vector<double> const&
TCHSNoStokesLocalAssembler<ShapeFunctionVelocity, ShapeFunctionPressure,
                           IntegrationMethod, VelocityDim>::
    getIntPtReactionRate(const double /*t*/,
                         GlobalVector const& /*current_solution*/,
                         NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
                         std::vector<double>& cache) const
{
    cache.clear();
    cache.reserve(_ip_data.size());

    auto const mat_id = _process_data.material_ids[_element.getID()];
    auto const& mat = _process_data.materials.at(mat_id);

    for (auto& ip_data : _ip_data)
    {
        cache.emplace_back(
            mat.reactive_solid->getOverallRate(*ip_data.reaction_rate));
    }
    return cache;
}

template <typename ShapeFunctionVelocity, typename ShapeFunctionPressure,
          typename IntegrationMethod, int VelocityDim>
std::vector<double> const&
TCHSNoStokesLocalAssembler<ShapeFunctionVelocity, ShapeFunctionPressure,
                           IntegrationMethod, VelocityDim>::
    getIntPtPorosity(const double /*t*/,
                     GlobalVector const& /*current_solution*/,
                     NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
                     std::vector<double>& cache) const
{
    cache.clear();
    cache.reserve(_ip_data.size());

    auto const mat_id = _process_data.material_ids[_element.getID()];
    auto const& mat = _process_data.materials.at(mat_id);

    for (auto& ip_data : _ip_data)
    {
        auto const r = interpolateXCoordinate<ShapeFunctionVelocity,
                                              ShapeMatricesTypeVelocity>(
            _element, ip_data.N_1);
        cache.emplace_back(mat.porosity->getPorosity(r));
    }
    return cache;
}

template <typename ShapeFunctionVelocity, typename ShapeFunctionPressure,
          typename IntegrationMethod, int VelocityDim>
void TCHSNoStokesLocalAssembler<
    ShapeFunctionVelocity, ShapeFunctionPressure, IntegrationMethod,
    VelocityDim>::preOutputConcrete(std::vector<double> const& local_x,
                                    const double t)
{
    auto const nodal_p = Block::mapVectorSegment(local_x, Block::P);
    auto const nodal_T = Block::mapVectorSegment(local_x, Block::T);
    auto const nodal_xmV = Block::mapVectorSegment(local_x, Block::X);

    SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    double cumul_hat_rho_SR = 0.0;
    double cumul_reaction_enthalpy = 0.0;
    double cumul_rho_SR = 0.0;
    double cumul_rho_GR = 0.0;
    Eigen::Matrix<double, VelocityDim, 1> cumul_lambda =
        Eigen::Matrix<double, VelocityDim, 1>::Zero();
    double cumul_cpS = 0.0;
    double cumul_cpG = 0.0;
    double cumul_volume = 0.0;

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto const& ip_data = _ip_data[ip];

        // shape functions and weights /////////////////////////////////////////
        x_position.setIntegrationPoint(ip);
        auto const& w = ip_data.integration_weight;
        auto const& N_1 = ip_data.N_1;

        auto const x_coord =
            interpolateXCoordinate<ShapeFunctionVelocity,
                                   ShapeMatricesTypeVelocity>(_element, N_1);

        // interpolate nodal values ////////////////////////////////////////////
        // TODO fix
        Eigen::Matrix<double, VelocityDim, 1> const v_Darcy =
            Eigen::Matrix<double, VelocityDim, 1>::Zero();
        double const p = N_1 * nodal_p;
        double const T = N_1 * nodal_T;
        double x_mV = N_1 * nodal_xmV;

        // material parameters /////////////////////////////////////////////////
        auto const mat_id = _process_data.material_ids[_element.getID()];
        auto const& mat = _process_data.materials.at(mat_id);

        double const M_R = mat.molar_mass_reactive;
        double const M_I = mat.molar_mass_inert;
        double p_V = p * Adsorption::AdsorptionReaction::getMolarFraction(
                             x_mV, M_R, M_I);

        // porosity
        auto const porosity = mat.porosity->getPorosity(x_coord);
        Eigen::Matrix<double, VelocityDim, 1> grad_porosity =
            Eigen::Matrix<double, VelocityDim, 1>::Zero();
        grad_porosity[0] = mat.porosity->getDPorosityDr(x_coord);

        // reaction
        if (mat.reaction_rate->computeReactionRate(
                _process_data.delta_t, p, T, p_V, *ip_data.reaction_rate,
                ip_data.reactive_solid_state.get(),
                ip_data.reaction_rate_data.get()))
        {
            x_mV = Adsorption::AdsorptionReaction::getMassFraction(p_V / p, M_R,
                                                                   M_I);
        }

        double const hat_rho_SR =
            mat.reactive_solid->getOverallRate(*ip_data.reaction_rate);
        auto const reaction_enthalpy = mat.reaction_rate->getHeatOfReaction(
            p_V, T, ip_data.reactive_solid_state.get());

        // fluid
        double const rho_GR = mat.fluid_density->getDensity(p, T, x_mV);
        double const c_pG = mat.fluid_heat_capacity->getHeatCapacity(T, x_mV);

        // fluid viscosity/friction
        double const mu = mat.fluid_viscosity->getViscosity(p, T, x_mV);
        double const Re_0 =
            mat.reynolds_number->getRe(t, rho_GR, v_Darcy.norm(), mu);

        // solid
        double const rho_SR =
            mat.reactive_solid->getSolidDensity(*ip_data.reactive_solid_state);
        double const c_pS =
            mat.solid_heat_capacity->getSpecificHeatCapacity(rho_SR, T);

        // fluid and solid
        auto const total_heat_conductivity =
            mat.heat_conductivity->getHeatConductivity(
                t, p, T, x_mV, x_coord, porosity, rho_GR, c_pG, Re_0,
                v_Darcy.norm(), _process_data.probed_velocity);

        cumul_hat_rho_SR += hat_rho_SR * w;
        cumul_reaction_enthalpy += reaction_enthalpy[0] * w;
        cumul_rho_SR += rho_SR * w;
        cumul_rho_GR += rho_GR * w;
        cumul_lambda += total_heat_conductivity.diagonal() * w;
        cumul_cpS += c_pS * w;
        cumul_cpG += c_pG * w;
        cumul_volume += w;
    }

    auto const id = _element.getID();
    (*_process_data.mesh_prop_cell_hat_rho_SR)[id] =
        cumul_hat_rho_SR / cumul_volume;
    (*_process_data.mesh_prop_cell_reaction_enthalpy)[id] =
        cumul_reaction_enthalpy / cumul_volume;
    (*_process_data.mesh_prop_cell_rho_SR)[id] = cumul_rho_SR / cumul_volume;
    (*_process_data.mesh_prop_cell_rho_GR)[id] = cumul_rho_GR / cumul_volume;
    for (int d = 0; d < VelocityDim; ++d)
    {
        (*_process_data.mesh_prop_cell_lambda)[VelocityDim * id + d] =
            cumul_lambda[d] / cumul_volume;
    }
    (*_process_data.mesh_prop_cell_cpS)[id] = cumul_cpS / cumul_volume;
    (*_process_data.mesh_prop_cell_cpG)[id] = cumul_cpG / cumul_volume;
}

}  // namespace TCHSNoStokes
}  // namespace ProcessLib
