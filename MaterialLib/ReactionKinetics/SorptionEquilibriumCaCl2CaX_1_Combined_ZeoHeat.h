/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#pragma once

#include "MaterialLib/ReactionKinetics/DubininPolanyi.h"
#include "SorptionEquilibriumCaCl2CaX_1_Combined.h"

namespace MaterialLib
{
class SorptionEquilibriumCaCl2CaX_1_Combined_ZeoHeat final
    : public SorptionEquilibriumCaCl2CaX_1_Combined
{
public:
    SorptionEquilibriumCaCl2CaX_1_Combined_ZeoHeat(const double rho_SR_dry,
                                                   const double salt_loading)
        : SorptionEquilibriumCaCl2CaX_1_Combined(rho_SR_dry, salt_loading)
    {
    }

    HeatOfReactionData getHeatOfReaction(
        const double p_Ads, const double T_Ads,
        ReactiveSolidState const* const state) const override;

private:
    MaterialLib::DubininPolanyi _caX_equilibrium{
        _rho_SR_dry,
        MaterialLib::createDubininPolanyiData("CaX80NoOutliers_Hauer")};
};

std::unique_ptr<SorptionEquilibriumCaCl2CaX_1_Combined_ZeoHeat>
createSorptionEquilibriumCaCl2CaX_1_Combined_ZeoHeat(
    BaseLib::ConfigTree const& config);

}  // namespace MaterialLib
