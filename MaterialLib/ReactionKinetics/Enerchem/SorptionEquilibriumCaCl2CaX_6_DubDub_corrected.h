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
class SorptionEquilibriumCaCl2CaX_6_DubDub_corrected final
    : public ReactionEquilibrium
{
public:
    SorptionEquilibriumCaCl2CaX_6_DubDub_corrected(
        const double rho_SR_dry, const double salt_loading,
        const double custom_enthalpy_factor);

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
    const double _custom_enthalpy_factor;
};

std::unique_ptr<SorptionEquilibriumCaCl2CaX_6_DubDub_corrected>
createSorptionEquilibriumCaCl2CaX_6_DubDub_corrected(
    BaseLib::ConfigTree const& config);

}  // namespace MaterialLib
