/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATERIALLIB_REACTIONRATESIMPLE_H
#define MATERIALLIB_REACTIONRATESIMPLE_H

#include "ReactionKinetics.h"
#include "ReactionRate.h"

namespace MaterialLib
{
class ReactionRateSimple final : public ReactionRate
{
public:
    ReactionRateSimple(double eq_react_lim,
                       double upper_pV_lim,
                       std::unique_ptr<ReactionKinetics>&& react_kin);

    bool isStateCompatible(ReactiveSolidState& state) const override;

    void preIteration(const unsigned iter) override;

    bool computeReactionRate(double const delta_t, double const p,
                             double const T, double& p_V, ReactiveSolidRate& qR,
                             ReactiveSolidState* solid_state,
                             ReactionRateData* reaction_rate_data) override;

    HeatOfReactionData getHeatOfReaction(
        const double p_V, const double T,
        const ReactiveSolidState* const state) const override;

    std::vector<NumLib::NamedFunction> getNamedFunctions() const override
    {
        return _react_kin->getNamedFunctions();
    }

    double const _equilibrium_reaction_vapour_pressure_limit;
    double const _upper_vapour_pressure_limit;
    std::unique_ptr<ReactionKinetics> _react_kin;
    unsigned _iter;
};

std::unique_ptr<ReactionRate> createReactionRateSimple(
    BaseLib::ConfigTree const& config);

}  // namespace MaterialLib

#endif  // MATERIALLIB_REACTIONRATESIMPLE_H
