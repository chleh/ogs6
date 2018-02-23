/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATERIALLIB_REACTIONRATEDATAUPDATEMULTIPLE_H
#define MATERIALLIB_REACTIONRATEDATAUPDATEMULTIPLE_H

#include "ReactionRateData.h"

namespace MaterialLib
{
struct ReactionRateDataUpdateMultiple final : public ReactionRateData
{
public:
    void preTimestep(ReactiveSolidState const& solid_state,
                     ReactiveSolidRate const& reaction_rate) override;

    Eigen::VectorXd solid_density_prev_ts;
    bool damped_in_current_ts;
    Eigen::VectorXd reaction_rate_max;
};

}  // namespace MaterialLib

#endif  // MATERIALLIB_REACTIONRATEDATAUPDATEMULTIPLE_H
