/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ReactionRateDataUpdateMultiple.h"
#include <cassert>

namespace MaterialLib
{
void ReactionRateDataUpdateMultiple::preTimestep(
    ReactiveSolidState const& solid_state,
    ReactiveSolidRate const& reaction_rate)
{
    solid_density_prev_ts = solid_state.conversion();
    reaction_rate_max = reaction_rate;

    damped_in_current_ts = false;
}

}  // namespace MaterialLib
