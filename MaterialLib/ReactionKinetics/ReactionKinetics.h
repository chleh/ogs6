/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATERIALLIB_REACTIONKINETICS_H
#define MATERIALLIB_REACTIONKINETICS_H

#include <memory>

#include "NumLib/NamedFunctionProvider.h"

#include "ReactionHeatModel.h"
#include "ReactiveSolidRate.h"

namespace BaseLib
{
class ConfigTree;
}  // namespace BaseLib

namespace MaterialLib
{
/*! Reaction kinetics for a reaction of type
 *  A(s) + B(g) <-> C(s).
 *
 * The reaction progress is computed from the density of the solid.
 */
class ReactionKinetics : public NumLib::NamedFunctionProvider,
                         public ReactionHeatModel
{
public:
    virtual bool isStateCompatible(ReactiveSolidState& state) const = 0;
    virtual ReactiveSolidRate getReactionRate(
        const double p_V, const double p, const double T,
        ReactiveSolidState const& solid_state) = 0;
};

std::unique_ptr<ReactionKinetics> createReactionKinetics(
    BaseLib::ConfigTree const& config);

}  // namespace MaterialLib

#endif  // MATERIALLIB_REACTIONKINETICS_H
