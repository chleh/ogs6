/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATERIALLIB_DUBININPOLANYI_H
#define MATERIALLIB_DUBININPOLANYI_H

#include "DubininPolanyiData.h"
#include "ReactionEquilibrium.h"

namespace MaterialLib
{
double waterEnthalpyOfEvaporation(const double T);

class DubininPolanyi : public ReactionEquilibrium
{
public:
    DubininPolanyi(const double rho_SR_dry,
                   DubininPolanyiData&& equilibrium_data)
        : _rho_SR_dry(rho_SR_dry), _dp_data(std::move(equilibrium_data))
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
    DubininPolanyiData _dp_data;
};

std::unique_ptr<DubininPolanyi> createDubininPolanyi(
    BaseLib::ConfigTree const& config);

}  // namespace MaterialLib

#endif  // MATERIALLIB_DUBININPOLANYI_H
