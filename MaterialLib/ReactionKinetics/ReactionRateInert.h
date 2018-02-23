/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ReactionRate.h"

namespace MaterialLib
{
class ReactionRateInert final : public ReactionRate
{
public:
    bool computeReactionRate(double const /*delta_t*/, double const /*p*/,
                             double const /*T*/, double& /*p_V*/,
                             ReactiveSolidRate& qR,
                             ReactiveSolidState* /*solid_state*/,
                             ReactionRateData* /*reaction_rate_data*/) override
    {
        qR.setZero();
        return false;
    }

    bool isStateCompatible(ReactiveSolidState& /*state*/) const override
    {
        return true;
    }

    HeatOfReactionData getHeatOfReaction(
        const double /*p_V*/, const double /*T*/,
        const ReactiveSolidState* const state) const override
    {
        return Eigen::VectorXd::Zero(state->conversion().size());
    }
};

}  // MaterialLib
