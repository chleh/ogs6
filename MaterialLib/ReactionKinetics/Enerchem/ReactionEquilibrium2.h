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

#include "NumLib/NamedFunctionProvider.h"

#include "MaterialLib/ReactionKinetics/HeatOfReactionData.h"
#include "MaterialLib/ReactionKinetics/ReactiveSolidState.h"

namespace BaseLib
{
class ConfigTree;
}  // namespace BaseLib

namespace MaterialLib
{
class ReactionEquilibrium2 : public NumLib::NamedFunctionProvider
{
public:
    virtual void getEquilibrium(const double p, const double T,
                                ReactiveSolidState const& solid_state,
                                ReactiveSolidState& equilibrium) const = 0;

    virtual HeatOfReactionData getHeatOfReaction(
        const double p_Ads, const double T_Ads,
        ReactiveSolidState const* const state) const = 0;

    virtual ~ReactionEquilibrium2() = default;
};

// indicates that equilibrium is based on zeolite and salt loading.
class ReactionEquilibrium2_LoadingLoading : public ReactionEquilibrium2
{
};

std::unique_ptr<ReactionEquilibrium2> createReactionEquilibrium2(
    BaseLib::ConfigTree const& config);

}  // namespace MaterialLib
