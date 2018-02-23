/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATERIALLIB_REACTIONRATEUPDATEONCE_H
#define MATERIALLIB_REACTIONRATEUPDATEONCE_H

#include "ReactionKinetics.h"
#include "ReactionRate.h"

namespace MaterialLib
{
class ReactionRateUpdateOnce final : public ReactionRate
{
public:
    ReactionRateUpdateOnce(double eq_react_lim, double hard_pV_lim,
                           double upper_pV_lim,
                           std::unique_ptr<ReactionKinetics>&& react_kin)
        : _equilibrium_reaction_vapour_pressure_limit(eq_react_lim),
          _hard_lower_vapour_pressure_limit(hard_pV_lim),
          _upper_vapour_pressure_limit(upper_pV_lim),
          _react_kin(std::move(react_kin))
    {
    }
    bool computeReactionRate(double const delta_t, double const p,
                             double const T, double& p_V, ReactiveSolidRate& qR,
                             ReactiveSolidState* solid_state,
                             ReactionRateData* reaction_rate_data) override;

    HeatOfReactionData getHeatOfReaction(
        const double p_V, const double T,
        ReactiveSolidState const* const state) const override
    {
        return _react_kin->getHeatOfReaction(p_V, T, state);
    }

    std::vector<NumLib::NamedFunction> getNamedFunctions() const override
    {
        return _react_kin->getNamedFunctions();
    }

    bool isStateCompatible(ReactiveSolidState& state) const override;

    void preIteration(const unsigned iter) override;
    void postIteration() override;

    std::unique_ptr<ReactionRateData> createReactionRateData() const override;

private:
    double const _equilibrium_reaction_vapour_pressure_limit;
    double const _hard_lower_vapour_pressure_limit;
    double const _upper_vapour_pressure_limit;
    std::unique_ptr<ReactionKinetics> _react_kin;
    unsigned _iter;
    bool _damped = false;
    // bool _damped_previous_iteration = false;
    unsigned _last_damped_iter = 0;
};

std::unique_ptr<ReactionRateUpdateOnce> createReactionRateUpdateOnce(
    BaseLib::ConfigTree const& config);

}  // namespace MaterialLib

#endif  // MATERIALLIB_REACTIONRATEUPDATEONCE_H
