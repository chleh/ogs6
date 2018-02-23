/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ReactionKinetics.h"
#include "ReactionRate.h"

namespace MaterialLib
{
class ReactionRateRaw final : public ReactionRate
{
public:
    ReactionRateRaw(std::unique_ptr<ReactionKinetics>&& react_kin);

    bool isStateCompatible(ReactiveSolidState& state) const override;

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

    std::unique_ptr<ReactionKinetics> _react_kin;
};

std::unique_ptr<ReactionRate> createReactionRateRaw(
    BaseLib::ConfigTree const& config);

}  // namespace MaterialLib
