/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATERIALLIB_REACTIONRATEDATA_H
#define MATERIALLIB_REACTIONRATEDATA_H

#include "ReactiveSolidRate.h"
#include "ReactiveSolidState.h"

namespace MaterialLib
{
class ReactionRateData
{
public:
    virtual void preTimestep(ReactiveSolidState const& solid_state,
                             ReactiveSolidRate const& reaction_rate) = 0;

    virtual ~ReactionRateData() = default;
};

}  // namespace MaterialLib

#endif  // MATERIALLIB_REACTIONRATEDATA_H
