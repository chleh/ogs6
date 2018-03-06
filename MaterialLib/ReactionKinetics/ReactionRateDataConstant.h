/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ReactionRateData.h"

namespace MaterialLib
{
struct ReactionRateDataConstant final : public ReactionRateData
{
public:
    void preTimestep(ReactiveSolidState const& solid_state,
                     ReactiveSolidRate const& reaction_rate) override;

    Eigen::VectorXd solid_density_prev_ts;
};

}  // namespace MaterialLib
