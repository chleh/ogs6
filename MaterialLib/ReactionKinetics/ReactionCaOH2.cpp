/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ReactionCaOH2.h"

#include <cassert>
#include <cmath>

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "MaterialLib/Adsorption/Adsorption.h"
#include "MaterialLib/PhysicalConstant.h"

// #define SIMPLE_KINETICS

namespace MaterialLib
{
const double ReactionCaOH2::_reaction_enthalpy = -1.12e+05;
const double ReactionCaOH2::_reaction_entropy = -143.5;
const double ReactionCaOH2::_M_react =
    MaterialLib::PhysicalConstant::MolarMass::Water;

const double ReactionCaOH2::_tol_l = 1e-4;
const double ReactionCaOH2::_tol_u = 1.0 - 1e-4;
const double ReactionCaOH2::_tol_rho = 0.1;

const double ReactionCaOH2::rho_low = 1665.0;
const double ReactionCaOH2::rho_up = 2200.0;

ReactiveSolidRate ReactionCaOH2::getReactionRate(
    const double p_V,
    const double /*p*/,
    const double T,
    ReactiveSolidState const& solid_state)
{
    assert(solid_state.conversion().size() == 1);
    auto const rho_SR = solid_state.conversion()[0];

    const double p_V_bar =
        std::max(p_V / 1e5, 1.e-3);  // convert Pa to bar; avoid illdefined log

    const double R = MaterialLib::PhysicalConstant::IdealGasConstant;

    double X_D =
        (rho_SR - rho_up - _tol_rho) / (rho_low - rho_up - 2.0 * _tol_rho);
    X_D = (X_D < 0.5)
              ? std::max(_tol_l, X_D)
              : std::min(X_D, _tol_u);  // constrain to interval [tol_l;tol_u]

    const double X_H = 1.0 - X_D;

    // calculate equilibrium
    // using the p_eq to calculate the T_eq - Clausius-Clapeyron
    const double T_eq =
        (_reaction_enthalpy / R) /
        ((_reaction_entropy / R) + std::log(p_V_bar));  // unit of p in bar
    // Alternative: Use T_s as T_eq and calculate p_eq - for Schaube kinetics
    const double p_eq =
        std::exp((_reaction_enthalpy / R) / T - (_reaction_entropy / R));

    // CaHydration():
    double dXdt;
// step 3, calculate dX/dt
#ifdef SIMPLE_KINETICS
    (void)p_eq;
    if (T < T_eq)  // hydration - simple model
#else
    if (p_V_bar > p_eq)  // hydration - Schaube model
#endif
    {
// X_H = max(tol_l,X_H); //lower tolerance to avoid oscillations at onset of
// hydration reaction. Set here so that no residual reaction rate occurs at end
// of hydration.
#ifdef SIMPLE_KINETICS               // this is from P. Schmidt
        const double x_react = 1.0;  // wrong: x_react = x_mV !!
        dXdt = -1.0 * (1.0 - X_H) * (T - T_eq) / T_eq * 0.2 * x_react;
#else  // this is from Schaube
        if (X_H == _tol_u || rho_SR == rho_up)
            dXdt = 0.0;
        else if ((T_eq - T) >= 50.0)
            dXdt = 13945.0 * std::exp(-89486.0 / R / T) *
                   std::pow(p_V_bar / p_eq - 1.0, 0.83) * 3.0 *
                   (X_D)*std::pow(-1.0 * log(X_D), 0.666);
        else
            dXdt = 1.0004e-34 * std::exp(5.3332e4 / T) *
                   std::pow(p_V_bar, 6.0) * (X_D);
#endif
    }
    else  // dehydration
    {
// X_D = max(tol_l,X_D); //lower tolerance to avoid oscillations at onset of
// dehydration reaction. Set here so that no residual reaction rate occurs at
// end of dehydration.
#ifdef SIMPLE_KINETICS  // this is from P. Schmidt
        dXdt = -1.0 * (1.0 - X_D) * (T - T_eq) / T_eq * 0.05;
#else
        if (X_D == _tol_u || rho_SR == rho_low)
            dXdt = 0.0;
        else if (X_D < 0.2)
            dXdt = -1.9425e12 * std::exp(-1.8788e5 / R / T) *
                   std::pow(1.0 - p_V_bar / p_eq, 3.0) * (X_H);
        else
            dXdt = -8.9588e9 * std::exp(-1.6262e5 / R / T) *
                   std::pow(1.0 - p_V_bar / p_eq, 3.0) * 2.0 *
                   std::pow(X_H, 0.5);
#endif
    }

    return Eigen::VectorXd::Constant(1, (rho_up - rho_low) * dXdt);
}

HeatOfReactionData ReactionCaOH2::getHeatOfReaction(
    const double /*p*/, const double /*T*/,
    const ReactiveSolidState* const /*state*/) const
{
    return (Eigen::VectorXd{1} << -_reaction_enthalpy / _M_react).finished();
}

// TODO use somewhere
#if 0
TESFEMReactionAdaptorCaOH2::TESFEMReactionAdaptorCaOH2(
    TESLocalAssemblerData const& data)
    : _d(data),
      _react(dynamic_cast<Adsorption::ReactionCaOH2&>(*data.ap.react_sys.get()))
{
    _ode_solver = MathLib::ODE::createODESolver<1>(_react.getOdeSolverConfig());
    // TODO invalidate config

    _ode_solver->setTolerance(1e-10, 1e-10);

    auto f = [this](const double /*t*/,
                    MathLib::ODE::MappedConstVector<1> const y,
                    MathLib::ODE::MappedVector<1>
                        ydot) -> bool {
        ydot[0] = _react.getReactionRate(y[0]);
        return true;
    };

    _ode_solver->setFunction(f, nullptr);
}

ReactionRate TESFEMReactionAdaptorCaOH2::initReaction(const unsigned int int_pt)
{
    // TODO if the first holds, the second also has to hold
    if (_d.ap.iteration_in_current_timestep > 1 ||
        _d.ap.number_of_try_of_iteration > 1)
    {
        return {_d.reaction_rate[int_pt], _d.solid_density[int_pt]};
    }

    // TODO: double check!
    // const double xv_NR  = SolidProp->non_reactive_solid_volume_fraction;
    // const double rho_NR = SolidProp->non_reactive_solid_density;
    const double xv_NR = 0.0;
    const double rho_NR = 0.0;

    const double t0 = 0.0;
    const double y0 =
        (_d.solid_density_prev_ts[int_pt] - xv_NR * rho_NR) / (1.0 - xv_NR);

    const double t_end = _d.ap.delta_t;

    _react.updateParam(_d.T, _d.p, _d.vapour_mass_fraction,
                       _d.solid_density_prev_ts[int_pt]);

    _ode_solver->setIC(t0, {y0});
    _ode_solver->preSolve();
    _ode_solver->solve(t_end);

    const double time_reached = _ode_solver->getTime();
    (void)time_reached;
    assert(std::abs(t_end - time_reached) <
           std::numeric_limits<double>::epsilon());

    auto const& y_new = _ode_solver->getSolution();
    auto const& y_dot_new = _ode_solver->getYDot(t_end, y_new);

    double rho_react;

    // cut off when limits are reached
    if (y_new[0] < _react.rho_low)
        rho_react = _react.rho_low;
    else if (y_new[0] > _react.rho_up)
        rho_react = _react.rho_up;
    else
        rho_react = y_new[0];

    return {y_dot_new[0] * (1.0 - xv_NR),
            (1.0 - xv_NR) * rho_react + xv_NR * rho_NR};
}
#endif

std::unique_ptr<ReactionCaOH2> createReactionCaOH2(
    BaseLib::ConfigTree const& config)
{
    config.checkConfigParameter("type", "CaOH2_Schaube");

    /*
    auto const k_rate = config.getConfigParameter<double>("k_rate");
    auto equil =
        createReactionEquilibrium(config.getConfigSubtree("equilibrium"));
        */

    return std::unique_ptr<ReactionCaOH2>(new ReactionCaOH2);
}

}  // namespace MaterialLib
