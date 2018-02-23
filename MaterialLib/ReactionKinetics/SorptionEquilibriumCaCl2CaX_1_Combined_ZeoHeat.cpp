/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "SorptionEquilibriumCaCl2CaX_1_Combined_ZeoHeat.h"
#include "BaseLib/ConfigTree.h"

namespace MaterialLib
{
HeatOfReactionData
SorptionEquilibriumCaCl2CaX_1_Combined_ZeoHeat::getHeatOfReaction(
    const double p_Ads,
    const double T_Ads,
    const ReactiveSolidState* const state) const
{
    return _caX_equilibrium.getHeatOfReaction(p_Ads, T_Ads, state);
}

std::unique_ptr<SorptionEquilibriumCaCl2CaX_1_Combined_ZeoHeat>
createSorptionEquilibriumCaCl2CaX_1_Combined_ZeoHeat(
    BaseLib::ConfigTree const& config)
{
    config.checkConfigParameter(
        "type", "SorptionEquilibriumCaCl2CaX_1_Combined_ZeoHeat");

    auto const rho_SR_dry =
        config.getConfigParameter<double>("sorbent_density_dry");
    auto const salt_loading =
        config.getConfigParameter<double>("salt_loading");  // mol/mol

    return std::unique_ptr<SorptionEquilibriumCaCl2CaX_1_Combined_ZeoHeat>(
        new SorptionEquilibriumCaCl2CaX_1_Combined_ZeoHeat(rho_SR_dry,
                                                           salt_loading));
}

}  // namespace MaterialLib
