/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "SorptionEquilibriumCaCl2CaX_5_DubDub.h"
#include "BaseLib/ConfigTree.h"
#include "BaseLib/Functional.h"

namespace MaterialLib
{
SorptionEquilibriumCaCl2CaX_5_DubDub::SorptionEquilibriumCaCl2CaX_5_DubDub(
    const double rho_SR_dry, const double salt_loading)
    : _alpha(1.0 - salt_loading / 15.0),
      _eq_00_salt{rho_SR_dry,
                  createDubininPolanyiData("CaX80NoOutliers_Hauer")},
      _eq_15_salt{rho_SR_dry, createDubininPolanyiData("15CaCl2CaX_Hauer")}
{
}

double SorptionEquilibriumCaCl2CaX_5_DubDub::getEquilibriumDensity(
    const double p_Ads, const double T_Ads)
{
    return _alpha * _eq_00_salt.getEquilibriumDensity(p_Ads, T_Ads) +
           (1.0 - _alpha) * _eq_15_salt.getEquilibriumDensity(p_Ads, T_Ads);
}

double SorptionEquilibriumCaCl2CaX_5_DubDub::getDEquilibriumDensityDp(
    const double p_Ads, const double T_Ads)
{
    return _alpha * _eq_00_salt.getDEquilibriumDensityDp(p_Ads, T_Ads) +
           (1.0 - _alpha) * _eq_15_salt.getDEquilibriumDensityDp(p_Ads, T_Ads);
}

HeatOfReactionData SorptionEquilibriumCaCl2CaX_5_DubDub::getHeatOfReaction(
    const double /*p_Ads*/, const double /*T_Ads*/,
    ReactiveSolidState const* const /*state*/) const
{
    // Approximation because the formula below yields too high values in order
    // to make the call A_of_W succeed internally.
    return (HeatOfReactionData(1) << 3.6e6).finished();
    /*
    return _alpha * _eq_00_salt.getHeatOfReaction(p_Ads, T_Ads, state) +
           (1.0 - _alpha) * _eq_15_salt.getHeatOfReaction(p_Ads, T_Ads, state);
           */
}

double SorptionEquilibriumCaCl2CaX_5_DubDub::getHeatOfReactionImpl(
    const double /*T_Ads*/, const double /*loading*/) const
{
    // Approximation because the formula below yields too high values in order
    // to make the call A_of_W succeed internally.
    return 3.6e6;
    /*
    return _alpha * _eq_00_salt.getHeatOfReactionImpl(T_Ads, loading) +
           (1.0 - _alpha) * _eq_15_salt.getHeatOfReactionImpl(T_Ads, loading);
           */
}

double SorptionEquilibriumCaCl2CaX_5_DubDub::getEquilibriumLoading(
    const double p_Ads, const double T_Ads) const
{
    return _alpha * _eq_00_salt.getEquilibriumLoading(p_Ads, T_Ads) +
           (1.0 - _alpha) * _eq_15_salt.getEquilibriumLoading(p_Ads, T_Ads);
}

double SorptionEquilibriumCaCl2CaX_5_DubDub::getLoading(
    const double adsorbent_density) const
{
    // Only need to call one of them, because both equilibrium models should
    // give the same results.
    return _eq_00_salt.getLoading(adsorbent_density);
}

std::vector<NumLib::NamedFunction>
SorptionEquilibriumCaCl2CaX_5_DubDub::getNamedFunctions() const
{
    return {{"equilibrium_loading",
             {"adsorptive_partial_pressure", "adsorbent_temperature"},
             BaseLib::easyBind(
                 &SorptionEquilibriumCaCl2CaX_5_DubDub::getEquilibriumLoading,
                 *this)},
            {"loading",
             {"adsorbent_density"},
             BaseLib::easyBind(
                 &SorptionEquilibriumCaCl2CaX_5_DubDub::getLoading, *this)},
            {"heat_of_reaction",
             {"adsorbent_temperature", "adsorbent_density"},
             BaseLib::easyBind(
                 &SorptionEquilibriumCaCl2CaX_5_DubDub::getHeatOfReactionImpl,
                 *this)}};
}

std::unique_ptr<SorptionEquilibriumCaCl2CaX_5_DubDub>
createSorptionEquilibriumCaCl2CaX_5_DubDub(BaseLib::ConfigTree const& config)
{
    config.checkConfigParameter("type", "SorptionEquilibriumCaCl2CaX_5_DubDub");

    auto const rho_SR_dry =
        config.getConfigParameter<double>("adsorbent_density_dry");
    auto const salt_loading = config.getConfigParameter<double>("salt_loading");
    return std::unique_ptr<SorptionEquilibriumCaCl2CaX_5_DubDub>(
        new SorptionEquilibriumCaCl2CaX_5_DubDub(rho_SR_dry, salt_loading));
}

}  // namespace MaterialLib
