/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ReactionRateRaw.h"

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "MathLib/Nonlinear/Root1D.h"

namespace MaterialLib
{
ReactionRateRaw::ReactionRateRaw(std::unique_ptr<ReactionKinetics>&& react_kin)
    : _react_kin(std::move(react_kin))
{
}

bool ReactionRateRaw::computeReactionRate(
    const double delta_t, double const p, const double T, double& p_V,
    ReactiveSolidRate& qR, ReactiveSolidState* solid_state,
    ReactionRateData* /*reaction_rate_data*/)
{
    qR = _react_kin->getReactionRate(p_V, p, T, *solid_state);
    *solid_state += delta_t * qR;

    return false;
}

HeatOfReactionData ReactionRateRaw::getHeatOfReaction(
    const double p_V, const double T,
    const ReactiveSolidState* const state) const
{
    return _react_kin->getHeatOfReaction(p_V, T, state);
}

bool ReactionRateRaw::isStateCompatible(ReactiveSolidState& state) const
{
    return _react_kin->isStateCompatible(state);
    // return dynamic_cast<ReactiveSolidStateOneComponent*>(&state) != nullptr;
}

std::unique_ptr<ReactionRate> createReactionRateRaw(
    BaseLib::ConfigTree const& config)
{
    config.checkConfigParameter("type", "Raw");

    auto react_kin =
        createReactionKinetics(config.getConfigSubtree("reaction_kinetics"));

    return std::unique_ptr<ReactionRateRaw>(
        new ReactionRateRaw(std::move(react_kin)));
}

}  // namespace MaterialLib
