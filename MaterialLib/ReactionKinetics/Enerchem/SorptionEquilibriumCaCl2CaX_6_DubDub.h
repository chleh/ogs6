/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MaterialLib/ReactionKinetics/ReactionEquilibrium.h"

namespace MaterialLib
{
class SorptionEquilibriumCaCl2CaX_6_DubDub : public ReactionEquilibrium
{
public:
    SorptionEquilibriumCaCl2CaX_6_DubDub(const double rho_SR_dry,
                                         const double salt_loading)
        : _rho_SR_dry(rho_SR_dry), _alpha(1.0 - salt_loading / 15.0)
    {
    }
    double getEquilibriumDensity(const double p_Ads,
                                 const double T_Ads) override;

    double getDEquilibriumDensityDp(const double p, const double T) override;

    HeatOfReactionData getHeatOfReaction(
        const double p_Ads, const double T_Ads,
        ReactiveSolidState const* const state) const override;

    std::vector<NumLib::NamedFunction> getNamedFunctions() const override;

    double getEntropy(const double T_Ads, const double A) const;

    virtual double getEquilibriumLoading(const double p_Ads,
                                         const double T_Ads) const;

    double getMaximumLoading(const double T) const;

    double getLoading(const double adsorbent_density) const;

    double getHeatOfReactionImpl(const double T_Ads,
                                 const double loading) const;

private:
    const double _rho_SR_dry;
    const double _alpha;
};

std::unique_ptr<SorptionEquilibriumCaCl2CaX_6_DubDub>
createSorptionEquilibriumCaCl2CaX_6_DubDub(BaseLib::ConfigTree const& config);

}  // namespace MaterialLib
