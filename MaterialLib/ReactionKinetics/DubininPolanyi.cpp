/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "DubininPolanyi.h"

#include <cmath>

#include <boost/math/tools/roots.hpp>

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "BaseLib/Functional.h"
#include "MaterialLib/PhysicalConstant.h"
#include "ReactionEquilibrium.h"

double A_of_W(const double W, MaterialLib::DubininPolanyiData const& dp_data)
{
    auto const A_min = 0.0;
    // TODO use a more general upper limit.
    auto const A_max =
        3500.0;  // the root of the CC for Na-X 13XBFK is already at ~3250 J/g

    auto const A_init = 0.5 * (A_max + A_min);
    int const digits = 8;
    std::uintmax_t max_iter_init = 100;
    std::uintmax_t max_iter = max_iter_init;

    auto const A = boost::math::tools::newton_raphson_iterate(
        [W, &dp_data](double A) {
            return std::make_pair(W - dp_data.characteristicCurve(A),
                                  -dp_data.characteristicCurveDeriv(A));
        },
        A_init, A_min, A_max, digits, max_iter);

    if (max_iter >= max_iter_init)
    {
        // if (A < 1e-3 && W > W_of_A(A, alpha))
        //     return 100; // DEBUG!!!!

        OGS_FATAL(
            "Dubinin-Polanyi A(W): Newton-Raphson scheme did not converge."
            " A=%g, searched W=%g, W(A)=%g",
            A, W, dp_data.characteristicCurve(A));
    }
    // INFO("A_of_W: A = %g", A);
    return A;
}

// Saturation pressure for water used in Nunez
double getEquilibriumVapourPressure(const double T_Ads)
{
    // critical T and p
    const double Tc = 647.3;    // K
    const double pc = 221.2e5;  // Pa
    // dimensionless T
    const double Tr = T_Ads / Tc;
    const double theta = 1. - Tr;
    // empirical constants
    const double c[] = {-7.69123, -26.08023, -168.17065, 64.23285, -118.96462,
                        4.16717,  20.97506,  1.0e9,      6.0};
    const double K[] = {c[0] * theta + c[1] * pow(theta, 2) +
                            c[2] * pow(theta, 3) + c[3] * pow(theta, 4) +
                            c[4] * pow(theta, 5),
                        1. + c[5] * theta + c[6] * pow(theta, 2)};

    const double exponent =
        K[0] / (K[1] * Tr) - theta / (c[7] * pow(theta, 2) + c[8]);
    return pc * exp(exponent);  // in Pa
}

double getPotential(const double p_Ads, const double T_Ads, const double M_Ads)
{
    const double p_S = getEquilibriumVapourPressure(T_Ads);
    const double A = MaterialLib::PhysicalConstant::IdealGasConstant * T_Ads *
                     std::log(p_S / p_Ads) / (M_Ads * 1.e3);  // in kJ/kg = J/g
    return A;
}

double getDPotentialDp(const double p_Ads, const double T_Ads,
                       const double M_Ads)
{
    return -MaterialLib::PhysicalConstant::IdealGasConstant * T_Ads /
           (M_Ads * 1.e3 * p_Ads);  // J / g / Pa
}

namespace MaterialLib
{
// Evaporation enthalpy of water from Nunez
double waterEnthalpyOfEvaporation(double T_Ads)  // in kJ/kg
{
    T_Ads -= 273.15;
    if (T_Ads <= 10.)
    {
        const double c[] = {2.50052e3,   -2.1068,     -3.57500e-1,
                            1.905843e-1, -5.11041e-2, 7.52511e-3,
                            -6.14313e-4, 2.59674e-5,  -4.421e-7};
        double hv = 0.;
        for (size_t i = 0; i < sizeof(c) / sizeof(c[0]); i++)
            hv += c[i] * pow(T_Ads, i);
        return hv;
    }
    else if (T_Ads <= 300.)
    {
        const double c[] = {2.50043e3,  -2.35209,    1.91685e-4,  -1.94824e-5,
                            2.89539e-7, -3.51199e-9, 2.06926e-11, -6.4067e-14,
                            8.518e-17,  1.558e-20,   -1.122e-22};
        double hv = 0.;
        for (size_t i = 0; i < sizeof(c) / sizeof(c[0]); i++)
            hv += c[i] * pow(T_Ads, i);
        return hv;
    }
    else
    {
        const double c[] = {2.99866e3, -3.1837e-3,  -1.566964e1,
                            -2.514e-6, 2.045933e-2, 1.0389e-8};
        return (
            (c[0] + c[2] * T_Ads + c[4] * pow(T_Ads, 2)) /
            (1. + c[1] * T_Ads + c[3] * pow(T_Ads, 2) + c[5] * pow(T_Ads, 3)));
    }
}

double DubininPolanyi::getEquilibriumDensity(const double p_Ads,
                                             const double T_Ads)
{
#if 0
    // TODO this is for testing.
    if (p_Ads <= 5.0)
        return _rho_SR_dry;
#endif

    const double C_eq = getEquilibriumLoading(p_Ads, T_Ads);
    return (C_eq + 1.0) * _rho_SR_dry;
}

double DubininPolanyi::getDEquilibriumDensityDp(const double p, const double T)
{
    auto const A = getPotential(p, T, _dp_data.M_Ads);            // J / g
    auto const dWdA = _dp_data.characteristicCurveDeriv(A);       // cm^3 / J
    auto const dAdp = getDPotentialDp(p, T, _dp_data.M_Ads);      // J / g / Pa
    return 1.e3 * _dp_data.getAdsorbateDensity(T) * dWdA * dAdp;  // 1 / Pa
}

double DubininPolanyi::getEquilibriumLoading(const double p_Ads,
                                             const double T_Ads) const
{
    if (p_Ads <= 0.0)
        return 0.0;  // TODO for testing/debugging

    double A = getPotential(p_Ads, T_Ads, _dp_data.M_Ads);
    if (A < 0.0)
        A = 0.0;

    const double C_eq =
        _dp_data.getAdsorbateDensity(T_Ads) * _dp_data.characteristicCurve(A);
    if (C_eq < 0.0)
        return 0.0;

    return C_eq;
}

double DubininPolanyi::getMaximumLoading(const double T) const
{
    const double C0 =
        _dp_data.getAdsorbateDensity(T) * _dp_data.characteristicCurve(0.0);
    return C0;
}

double DubininPolanyi::getLoading(const double adsorbent_density) const
{
    return adsorbent_density / _rho_SR_dry - 1.0;
}

HeatOfReactionData DubininPolanyi::getHeatOfReaction(
    const double /*p_Ads*/, const double T_Ads,
    const ReactiveSolidState* const state) const
{
    assert(state->conversion().size() == 1);
    return (Eigen::VectorXd{1}
            << getHeatOfReactionImpl(T_Ads, state->conversion()[0]))
        .finished();
}

double DubininPolanyi::getHeatOfReactionImpl(const double T_Ads,
                                             const double loading) const
{
    if (_dp_data.custom_enthalpy_method)
        return _dp_data.custom_enthalpy_method(T_Ads, loading);

    auto const loading2 = getLoading(loading);

    const double W = loading2 / _dp_data.getAdsorbateDensity(T_Ads);

    const double A = A_of_W(W, _dp_data);

    return (waterEnthalpyOfEvaporation(T_Ads) + A -
            T_Ads * getEntropy(T_Ads, A)) *
           1000.0;  // in J/kg
}

// Calculate sorption entropy
double DubininPolanyi::getEntropy(const double T_Ads, const double A) const
{
    const double W = _dp_data.characteristicCurve(A);
    const double dWdA = _dp_data.characteristicCurveDeriv(A);
    const double dAdlnW = W / dWdA;

    return dAdlnW * _dp_data.getAdsorbateThermalExpansion(T_Ads);
}

std::vector<NumLib::NamedFunction> DubininPolanyi::getNamedFunctions() const
{
    return {{"equilibrium_loading",
             {"adsorptive_partial_pressure", "adsorbent_temperature"},
             BaseLib::easyBind(&DubininPolanyi::getEquilibriumLoading, *this)},
            {"loading",
             {"adsorbent_density"},
             BaseLib::easyBind(&DubininPolanyi::getLoading, *this)},
            {"heat_of_reaction",
             {"adsorbent_temperature", "adsorbent_density"},
             BaseLib::easyBind(&DubininPolanyi::getHeatOfReactionImpl, *this)}};
}

std::unique_ptr<DubininPolanyi> createDubininPolanyi(
    BaseLib::ConfigTree const& config)
{
    config.checkConfigParameter("type", "DubininPolanyi");

    auto const rho_SR_dry =
        config.getConfigParameter<double>("adsorbent_density_dry");
    auto data =
        createDubininPolanyiData(config.getConfigSubtree("equilibrium_data"));
    return std::unique_ptr<DubininPolanyi>(
        new DubininPolanyi(rho_SR_dry, std::move(data)));
}

}  // namespace MaterialLib
