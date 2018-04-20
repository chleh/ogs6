/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "SorptionEquilibriumCaCl2CaX_6_DubDub_corrected.h"
#include <boost/math/tools/roots.hpp>
#include <cmath>
#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "BaseLib/Functional.h"
#include "MaterialLib/PhysicalConstant.h"
#include "MaterialLib/ReactionKinetics/DubininPolanyi.h"

namespace
{
static const MaterialLib::DubininPolanyiData data_00 =
    MaterialLib::createDubininPolanyiData("CaX80NoOutliers_Hauer");

static const MaterialLib::DubininPolanyiData data_15 =
    MaterialLib::createDubininPolanyiData("15CaCl2CaX_Hauer");

static double W_of_A(const double A, const double alpha)
{
    return alpha * data_00.characteristicCurve(A) +
           (1.0 - alpha) * data_15.characteristicCurve(A);
}

static double dW_dA(const double A, const double alpha)
{
    return alpha * data_00.characteristicCurveDeriv(A) +
           (1.0 - alpha) * data_15.characteristicCurveDeriv(A);
}

#if 0
double A_of_W(const double W, const double alpha)
{
    INFO("A_of_W: W = %g", W);
    double A = 1500.0;

    const int max_iter = 100;

    const double abs_tol = 1e-10;
    const double rel_tol = 1e-8;

    for (int i = 0; i < max_iter; ++i) {
        const double f_A = W - W_of_A(A, alpha);
        const double minus_fp_A = dW_dA(A, alpha);
        const double delta_A = f_A / minus_fp_A;
        A += delta_A;

        if (A < 0.0)
            A = 0.0;
        // else if (A > 2500.0)
        //     A = 2500.0;
        else if (std::abs(delta_A) < abs_tol || std::abs(delta_A / W) < rel_tol)
            return A;
    }

    OGS_FATAL(
        "Dubinin-Polanyi A(W): Newton-Raphson scheme did not converge."
        " A=%g, W=%g",
        A, W);
}
#else
double A_of_W(const double W, const double alpha)
{
    // INFO("A_of_W: W = %g", W);
    auto const A_min = 0.0;
    auto const A_max = 2548.0;
    // CC data_00 has its root at 2548.740080014731,
    // and CC data_15 has its root at 2700.941049348698

    auto const A_init = 0.5 * (A_max + A_min);
    int const digits = 8;
    std::uintmax_t max_iter_init = 100;
    std::uintmax_t max_iter = max_iter_init;

    auto const A = boost::math::tools::newton_raphson_iterate(
        [W, alpha](double A) {
            return std::make_pair(W - W_of_A(A, alpha), -dW_dA(A, alpha));
        },
        A_init, A_min, A_max, digits, max_iter);

    if (max_iter >= max_iter_init)
    {
        // if (A < 1e-3 && W > W_of_A(A, alpha))
        //     return 100; // DEBUG!!!!

        OGS_FATAL(
            "Dubinin-Polanyi A(W): Newton-Raphson scheme did not converge."
            " A=%g, searched W=%g, W(A)=%g",
            A, W, W_of_A(A, alpha));
    }
    // INFO("A_of_W: A = %g", A);
    return A;
}
#endif

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
}

double getSaltLoadingAlpha(const double salt_loading_mole_fraction)
{
    auto const M_CaX = 13273;
    auto const M_CaCl2 = 110.980;
    auto const salt_loading_mass_fraction =
        salt_loading_mole_fraction * M_CaCl2 / M_CaX;
    auto const salt_content_mass_fraction =
        salt_loading_mass_fraction / (1.0 + salt_loading_mass_fraction);
    return salt_content_mass_fraction;
}

namespace MaterialLib
{
SorptionEquilibriumCaCl2CaX_6_DubDub_corrected::
    SorptionEquilibriumCaCl2CaX_6_DubDub_corrected(const double rho_SR_dry,
                                                   const double salt_loading)
    : _rho_SR_dry(rho_SR_dry),
      _alpha(1.0 - getSaltLoadingAlpha(salt_loading) / 15.0)
{
}

double SorptionEquilibriumCaCl2CaX_6_DubDub_corrected::getEquilibriumDensity(
    const double p_Ads, const double T_Ads)
{
#if 0
    // TODO this is for testing.
    if (p_Ads <= 5.0)
        return _rho_SR_dry;
#endif

    const double C_eq = getEquilibriumLoading(p_Ads, T_Ads);
    return (C_eq + 1.0) * _rho_SR_dry;
}

double SorptionEquilibriumCaCl2CaX_6_DubDub_corrected::getDEquilibriumDensityDp(
    const double p, const double T)
{
    auto const A = getPotential(p, T, data_00.M_Ads);            // J / g
    auto const dWdA = dW_dA(A, _alpha);                          // cm^3 / J
    auto const dAdp = getDPotentialDp(p, T, data_00.M_Ads);      // J / g / Pa
    return 1.e3 * data_00.getAdsorbateDensity(T) * dWdA * dAdp;  // 1 / Pa
}

double SorptionEquilibriumCaCl2CaX_6_DubDub_corrected::getEquilibriumLoading(
    const double p_Ads, const double T_Ads) const
{
    if (p_Ads <= 0.0)
        return 0.0;  // TODO for testing/debugging

    double A = getPotential(p_Ads, T_Ads, data_00.M_Ads);
    if (A < 0.0)
        A = 0.0;

    const double C_eq = data_00.getAdsorbateDensity(T_Ads) * W_of_A(A, _alpha);
    if (C_eq < 0.0)
        return 0.0;

    return C_eq;
}

double SorptionEquilibriumCaCl2CaX_6_DubDub_corrected::getMaximumLoading(
    const double T) const
{
    const double C0 = data_00.getAdsorbateDensity(T) * W_of_A(0.0, _alpha);
    return C0;
}

double SorptionEquilibriumCaCl2CaX_6_DubDub_corrected::getLoading(
    const double adsorbent_density) const
{
    return adsorbent_density / _rho_SR_dry - 1.0;
}

HeatOfReactionData
SorptionEquilibriumCaCl2CaX_6_DubDub_corrected::getHeatOfReaction(
    const double /*p_Ads*/, const double T_Ads,
    const ReactiveSolidState* const state) const
{
    assert(state->conversion().size() == 1);
    return (Eigen::VectorXd{1}
            << getHeatOfReactionImpl(T_Ads, state->conversion()[0]))
        .finished();
}

double SorptionEquilibriumCaCl2CaX_6_DubDub_corrected::getHeatOfReactionImpl(
    const double T_Ads, const double loading) const
{
    auto const loading2 = getLoading(loading);

    const double W = loading2 / data_00.getAdsorbateDensity(T_Ads);

    const double A = A_of_W(W, _alpha);

    // the entropic term is omitted here, because it only contributes
    // on average < 4% and misbehaves for some loading values
    return (waterEnthalpyOfEvaporation(T_Ads) + A
            // - T_Ads * getEntropy(T_Ads, A)
            ) *
           1000.0;  // in J/kg
}

// Calculate sorption entropy
double SorptionEquilibriumCaCl2CaX_6_DubDub_corrected::getEntropy(
    const double T_Ads, const double A) const
{
    const double W = W_of_A(A, _alpha);
    const double dWdA = dW_dA(A, _alpha);
    const double dAdlnW = W / dWdA;

    return dAdlnW * data_00.getAdsorbateThermalExpansion(T_Ads);
}

std::vector<NumLib::NamedFunction>
SorptionEquilibriumCaCl2CaX_6_DubDub_corrected::getNamedFunctions() const
{
    return {{"equilibrium_loading",
             {"adsorptive_partial_pressure", "adsorbent_temperature"},
             BaseLib::easyBind(&SorptionEquilibriumCaCl2CaX_6_DubDub_corrected::
                                   getEquilibriumLoading,
                               *this)},
            {"loading",
             {"adsorbent_density"},
             BaseLib::easyBind(
                 &SorptionEquilibriumCaCl2CaX_6_DubDub_corrected::getLoading,
                 *this)},
            {"heat_of_reaction",
             {"adsorbent_temperature", "adsorbent_density"},
             BaseLib::easyBind(&SorptionEquilibriumCaCl2CaX_6_DubDub_corrected::
                                   getHeatOfReactionImpl,
                               *this)}};
}

std::unique_ptr<SorptionEquilibriumCaCl2CaX_6_DubDub_corrected>
createSorptionEquilibriumCaCl2CaX_6_DubDub_corrected(
    BaseLib::ConfigTree const& config)
{
    config.checkConfigParameter(
        "type", "SorptionEquilibriumCaCl2CaX_6_DubDub_corrected");

    auto const rho_SR_dry =
        config.getConfigParameter<double>("adsorbent_density_dry");
    auto const salt_loading = config.getConfigParameter<double>("salt_loading");

    return std::unique_ptr<SorptionEquilibriumCaCl2CaX_6_DubDub_corrected>(
        new SorptionEquilibriumCaCl2CaX_6_DubDub_corrected(rho_SR_dry,
                                                           salt_loading));
}

}  // namespace MaterialLib
