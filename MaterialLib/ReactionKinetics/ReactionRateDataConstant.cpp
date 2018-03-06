/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ReactionRateDataConstant.h"
#include <cassert>

namespace MaterialLib
{
void ReactionRateDataConstant::preTimestep(
    ReactiveSolidState const& solid_state,
    ReactiveSolidRate const& /*reaction_rate*/)
{
    solid_density_prev_ts = solid_state.conversion();
}

}  // namespace MaterialLib
