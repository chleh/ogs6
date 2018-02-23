/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MaterialLib/ReactionKinetics/ReactionKinetics.h"
#include "MaterialLib/ReactionKinetics/ReactionRate.h"

namespace MaterialLib
{
class ReactionRateUpdateMultipleEnerchem final : public ReactionRate
{
public:
    ReactionRateUpdateMultipleEnerchem(
        double eq_react_lim, double upper_pV_lim, unsigned max_updates,
        std::unique_ptr<ReactionKinetics>&& react_kin)
        : _equilibrium_reaction_vapour_pressure_limit(eq_react_lim),
          _upper_vapour_pressure_limit(upper_pV_lim),
          _maximum_reaction_rate_updates(max_updates),
          _react_kin(std::move(react_kin))
    {
    }
    bool computeReactionRate(double const delta_t, double const /*p*/,
                             double const T, double& p_V, ReactiveSolidRate& qR,
                             ReactiveSolidState* solid_state,
                             ReactionRateData* reaction_rate_data) override;

    HeatOfReactionData getHeatOfReaction(
        const double p_V, const double T,
        ReactiveSolidState const* const state) const override
    {
        return _react_kin->getHeatOfReaction(p_V, T, state);
    }

    std::vector<NumLib::NamedFunction> getNamedFunctions() const override;

    bool isStateCompatible(ReactiveSolidState& state) const override;

    void preIteration(const unsigned iter) override;
    void postIteration() override;

    std::unique_ptr<ReactionRateData> createReactionRateData() const override;

private:
    double computeEquilibriumVapourPressure(
        const double p, const double T,
        const ReactiveSolidState& solid_state) const;

    double const _equilibrium_reaction_vapour_pressure_limit;
    double const _upper_vapour_pressure_limit;
    unsigned const _maximum_reaction_rate_updates;

    std::unique_ptr<ReactionKinetics> _react_kin;

    unsigned _iter;
    bool _damped = false;
    // bool _damped_previous_iteration = false;
    unsigned _last_damped_iter = 0;
};

std::unique_ptr<ReactionRateUpdateMultipleEnerchem>
createReactionRateUpdateMultipleEnerchem(BaseLib::ConfigTree const& config);

}  // namespace MaterialLib
