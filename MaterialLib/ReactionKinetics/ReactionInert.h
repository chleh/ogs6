/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATERIALLIB_REACTIONINERT_H
#define MATERIALLIB_REACTIONINERT_H

#include "BaseLib/Error.h"

#include "ReactionKinetics.h"

namespace MaterialLib
{
class ReactionInert final : public ReactionKinetics
{
public:
    ReactiveSolidRate getReactionRate(
        const double /*p_V*/, const double /*p*/, const double /*T*/,
        ReactiveSolidState const& solid_state) override
    {
        return Eigen::VectorXd::Zero(solid_state.conversion().size());
    }

    HeatOfReactionData getHeatOfReaction(
        const double /*p_V*/, const double /*T*/,
        const ReactiveSolidState* const state) const
    {
        return Eigen::VectorXd::Zero(state->conversion().size());
        ;
    }

    bool isStateCompatible(ReactiveSolidState& /*state*/) const override
    {
        return true;
    }
};

}  // namespace MaterialLib
#endif  // MATERIALLIB_REACTIONINERT_H
