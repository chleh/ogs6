/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include "HeatOfReactionData.h"
#include "ReactiveSolidState.h"

namespace BaseLib
{
class ConfigTree;
}  // namespace BaseLib

namespace MaterialLib
{
class ReactionHeatModel
{
public:
    virtual HeatOfReactionData getHeatOfReaction(
        const double p_V, const double T,
        ReactiveSolidState const* const state) const = 0;

    ~ReactionHeatModel() = default;
};

}  // namespace MaterialLib
