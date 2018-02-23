/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

// #include "MaterialLib/ReactionKinetics/ReactionEquilibrium.h"
#include "MaterialLib/ReactionKinetics/DubininPolanyi.h"

namespace MaterialLib
{
class SorptionEquilibriumCaCl2CaX_5_DubDub : public ReactionEquilibrium
{
public:
    SorptionEquilibriumCaCl2CaX_5_DubDub(const double rho_SR_dry,
                                         const double salt_loading);

    double getEquilibriumDensity(const double p_Ads,
                                 const double T_Ads) override;

    double getDEquilibriumDensityDp(const double p, const double T) override;

    HeatOfReactionData getHeatOfReaction(
        const double p_Ads, const double T_Ads,
        ReactiveSolidState const* const state) const override;

    std::vector<NumLib::NamedFunction> getNamedFunctions() const override;

    double getEquilibriumLoading(const double p_Ads, const double T_Ads) const;

    double getLoading(const double adsorbent_density) const;

    double getHeatOfReactionImpl(const double T_Ads,
                                 const double loading) const;

protected:
    const double _alpha;  //!< Convex combination factor of the two equilibria
                          //! based on salt loading.
    DubininPolanyi _eq_00_salt;
    DubininPolanyi _eq_15_salt;
};

std::unique_ptr<SorptionEquilibriumCaCl2CaX_5_DubDub>
createSorptionEquilibriumCaCl2CaX_5_DubDub(BaseLib::ConfigTree const& config);

}  // namespace MaterialLib
