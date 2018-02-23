/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATERIALLIB_REACTIONRATEDATAUPDATEONCE_H
#define MATERIALLIB_REACTIONRATEDATAUPDATEONCE_H

#include "ReactionRateData.h"

namespace MaterialLib
{
struct ReactionRateDataUpdateOnce final : public ReactionRateData
{
public:
    void preTimestep(ReactiveSolidState const& solid_state,
                     ReactiveSolidRate const& reaction_rate) override;

    Eigen::VectorXd solid_density_prev_ts;
};

}  // namespace MaterialLib

#endif  // MATERIALLIB_REACTIONRATEDATAUPDATEONCE_H
