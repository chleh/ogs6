#include "ReactionMnO.h"

#include <cassert>
#include <cmath>

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "MaterialLib/Adsorption/Adsorption.h"
#include "MaterialLib/PhysicalConstant.h"

// #define SIMPLE_KINETICS

namespace MaterialLib
{
const double ReactionMnO::_reaction_enthalpy = -1.12e+05;
const double ReactionMnO::_reaction_entropy = -143.5;
const double ReactionMnO::_M_react =
    MaterialLib::PhysicalConstant::MolarMass::Water;

const double ReactionMnO::_tol_l = 1e-4;
const double ReactionMnO::_tol_u = 1.0 - 1e-4;
const double ReactionMnO::_tol_rho = 0.1;

const double ReactionMnO::rho_low = 1665.0;
const double ReactionMnO::rho_up = 2200.0;

ReactiveSolidRate ReactionMnO::getReactionRate(
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

HeatOfReactionData ReactionMnO::getHeatOfReaction(
    const double /*p*/, const double /*T*/,
    const ReactiveSolidState* const /*state*/) const
{
    return (Eigen::VectorXd{1} << -_reaction_enthalpy / _M_react).finished();
}

std::unique_ptr<ReactionMnO> createReactionMnO(
    BaseLib::ConfigTree const& config)
{
    config.checkConfigParameter("type", "MnO");

    /*
    auto const k_rate = config.getConfigParameter<double>("k_rate");
    auto equil =
        createReactionEquilibrium(config.getConfigSubtree("equilibrium"));
        */

    return std::unique_ptr<ReactionMnO>(new ReactionMnO);
}

}  // namespace MaterialLib
