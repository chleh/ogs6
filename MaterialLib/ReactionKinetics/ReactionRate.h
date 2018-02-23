/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATERIALLIB_REACTIONRATE_H
#define MATERIALLIB_REACTIONRATE_H

#include <memory>

#include "NumLib/NamedFunctionProvider.h"

#include "ReactionHeatModel.h"
#include "ReactionRateData.h"
#include "ReactiveSolidRate.h"

namespace BaseLib
{
class ConfigTree;
}  // namespace BaseLib

namespace MaterialLib
{
class ReactionRate : public NumLib::NamedFunctionProvider,
                     public ReactionHeatModel
{
public:
    virtual bool isStateCompatible(ReactiveSolidState& state) const = 0;
    virtual void preIteration(const unsigned iter) { (void)iter; }
    virtual void postIteration() {}
    virtual bool computeReactionRate(double const delta_t, double const p,
                                     double const T, double& p_V,
                                     ReactiveSolidRate& qR,
                                     ReactiveSolidState* solid_state,
                                     ReactionRateData* reaction_rate_data) = 0;

    virtual std::unique_ptr<ReactionRateData> createReactionRateData() const
    {
        return nullptr;
    }
};

std::unique_ptr<ReactionRate> createReactionRate(
    BaseLib::ConfigTree const& config);

}  // namespace MaterialLib

#endif  // MATERIALLIB_REACTIONRATE_H
