/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#pragma once

#include "ReactionEquilibrium.h"

namespace MaterialLib
{
class SorptionEquilibriumCaCl2CaX_1_Combined : public ReactionEquilibrium
{
public:
    SorptionEquilibriumCaCl2CaX_1_Combined(const double rho_SR_dry,
                                           const double salt_loading)
        : _rho_SR_dry(rho_SR_dry), _salt_loading(salt_loading)
    {
    }

    double getEquilibriumDensity(const double p_Ads,
                                 const double T_Ads) override;

    double getDEquilibriumDensityDp(const double p, const double T) override;

    HeatOfReactionData getHeatOfReaction(
        const double p_Ads, const double T_Ads,
        ReactiveSolidState const* const state) const override;

    std::vector<NumLib::NamedFunction> getNamedFunctions() const override;

    double getEquilibriumLoading(const double p_Ads, const double T_Ads) const;

protected:
    const double _rho_SR_dry;
    const double _salt_loading;  //!< in mol/mol
};

std::unique_ptr<SorptionEquilibriumCaCl2CaX_1_Combined>
createSorptionEquilibriumCaCl2CaX_1_Combined(BaseLib::ConfigTree const& config);

}  // namespace MaterialLib
