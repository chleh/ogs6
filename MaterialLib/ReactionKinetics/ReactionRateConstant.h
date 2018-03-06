/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "BaseLib/ConfigTree.h"
#include "ReactionRate.h"
#include "ReactionRateDataConstant.h"

namespace MaterialLib
{
class ReactionRateConstant final : public ReactionRate
{
public:
    ReactionRateConstant(double rate, double heat_of_reaction)
        : _rate(rate), _heat_of_reaction(heat_of_reaction)
    {
    }

    bool computeReactionRate(double const delta_t, double const /*p*/,
                             double const /*T*/, double& /*p_V*/,
                             ReactiveSolidRate& qR,
                             ReactiveSolidState* solid_state,
                             ReactionRateData* reaction_rate_data) override
    {
        qR.setConstant(_rate);

        assert(dynamic_cast<ReactionRateDataConstant*>(reaction_rate_data) !=
               nullptr);
        auto& d = *static_cast<ReactionRateDataConstant*>(reaction_rate_data);

        auto const& rhoSR_prev_ts = d.solid_density_prev_ts;
        solid_state->conversion() = rhoSR_prev_ts + qR * delta_t;
        return false;
    }

    bool isStateCompatible(ReactiveSolidState& /*state*/) const override
    {
        return true;
    }

    HeatOfReactionData getHeatOfReaction(
        const double /*p_V*/, const double /*T*/,
        const ReactiveSolidState* const state) const override
    {
        return Eigen::VectorXd::Constant(state->conversion().size(),
                                         _heat_of_reaction);
    }

    std::unique_ptr<ReactionRateData> createReactionRateData() const override
    {
        return std::make_unique<ReactionRateDataConstant>();
    }

private:
    double const _rate;
    double const _heat_of_reaction;
};

std::unique_ptr<ReactionRateConstant> createReactionRateConstant(
    BaseLib::ConfigTree const& config)
{
    config.checkConfigParameter("type", "Constant");
    auto const rate = config.getConfigParameter<double>("rate");
    auto const heat_of_reaction =
        config.getConfigParameter<double>("heat_of_reaction");
    return std::make_unique<ReactionRateConstant>(rate, heat_of_reaction);
}

}  // MaterialLib
