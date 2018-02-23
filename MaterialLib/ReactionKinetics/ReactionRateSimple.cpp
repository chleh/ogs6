/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ReactionRateSimple.h"

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "MathLib/Nonlinear/Root1D.h"

#include "ReactiveSolidStateOneComponent.h"

namespace MaterialLib
{
ReactionRateSimple::ReactionRateSimple(
    double eq_react_lim, double upper_pV_lim,
    std::unique_ptr<ReactionKinetics>&& react_kin)
    : _equilibrium_reaction_vapour_pressure_limit(eq_react_lim),
      _upper_vapour_pressure_limit(upper_pV_lim),
      _react_kin(std::move(react_kin))
{
}

bool ReactionRateSimple::computeReactionRate(
    const double delta_t, double const p, const double T, double& p_V,
    ReactiveSolidRate& qR, ReactiveSolidState* solid_state,
    ReactionRateData* /*reaction_rate_data*/)
{
    if (p_V < _equilibrium_reaction_vapour_pressure_limit)
    {
        // For low vapour mass fraction compute reaction equilibrium.
        auto const f = [&](double curr_p_V) -> double {
            auto qR = _react_kin->getReactionRate(curr_p_V, p, T, *solid_state);
            assert(qR.size() == 1);
            return qR[0];
        };
        // Find root of reaction rate (== reaction equilibrium).
        auto rf =
            MathLib::Nonlinear::makeRegulaFalsi<MathLib::Nonlinear::Pegasus>(
                f, 1e-6, _upper_vapour_pressure_limit);
        if (rf.step(5))
        {
            p_V = rf.getResult();
        }
        else
        {
            auto p_q = rf.getEndWithLowerValue();
            // TODO make eaction rate limit configurable
            if (std::abs(p_q.second) < 1e-4)
                p_V = p_q.first;
            else
                OGS_FATAL(
                    "Regula Falsi failed. The signs of function values at both "
                    "interval ends are the same.");
        }

        // Can be interpreted as average of qR at the beginning of the timestep
        // and 0.0 (i.e. equilibrium reached) at the end of the timestep.
        qR /= 2.0;
        *solid_state += delta_t * qR;

        return true;
    }
    else if (_iter == 1)
    {
        qR = _react_kin->getReactionRate(p_V, p, T, *solid_state);
        *solid_state += delta_t * qR;
    }

    return false;
}

HeatOfReactionData ReactionRateSimple::getHeatOfReaction(
    const double p_V, const double T,
    const ReactiveSolidState* const state) const
{
    return _react_kin->getHeatOfReaction(p_V, T, state);
}

bool ReactionRateSimple::isStateCompatible(ReactiveSolidState& state) const
{
    return dynamic_cast<ReactiveSolidStateOneComponent*>(&state) != nullptr;
}

void ReactionRateSimple::preIteration(const unsigned iter)
{
    _iter = iter;
}

std::unique_ptr<ReactionRate> createReactionRateSimple(
    BaseLib::ConfigTree const& config)
{
    config.checkConfigParameter("type", "Simple");

    auto const eq_react_lim = config.getConfigParameter<double>(
        "equilibrium_reaction_vapour_pressure_limit");

    auto react_kin =
        createReactionKinetics(config.getConfigSubtree("reaction_kinetics"));

    auto const upper_pV_lim =
        config.getConfigParameter<double>("upper_vapour_pressure_limit");

    return std::unique_ptr<ReactionRateSimple>(new ReactionRateSimple(
        eq_react_lim, upper_pV_lim, std::move(react_kin)));
}

}  // namespace MaterialLib
