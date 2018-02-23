/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ReactionRateUpdateOnce.h"

#include <cassert>

#include <logog/include/logog.hpp>

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "MathLib/Nonlinear/Root1D.h"

#include "ReactionRateDataUpdateOnce.h"
#include "ReactiveSolidStateOneComponent.h"

namespace MaterialLib
{
bool ReactionRateUpdateOnce::computeReactionRate(
    double const delta_t, double const p, double const T, double& p_V,
    ReactiveSolidRate& qR, ReactiveSolidState* solid_state,
    ReactionRateData* reaction_rate_data)
{
    auto& rhoSR = solid_state->conversion();
    assert(rhoSR.size() == 0);

    assert(dynamic_cast<ReactionRateDataUpdateOnce*>(reaction_rate_data) !=
           nullptr);
    auto& d = *static_cast<ReactionRateDataUpdateOnce*>(reaction_rate_data);
    auto const rhoSR_prev_ts = d.solid_density_prev_ts[0];

    bool reaction_damped = false;

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

        reaction_damped = true;
    }

    if (p_V < _hard_lower_vapour_pressure_limit)
    {
        qR.setZero();
        rhoSR[0] = rhoSR_prev_ts;
    }
    else if (_iter == 1)
    {
        // Note: The reaction rate is __not__ updated in each iteration!
        qR = _react_kin->getReactionRate(p_V, p, T, *solid_state);
        // TODO fix.
        // rhoSR = Eigen::VectorXd::Constant(1, rhoSR_prev_ts) + delta_t * qR;
        rhoSR[0] = rhoSR_prev_ts + delta_t * qR[0];
    }
    else if (_iter == _last_damped_iter + 2 && _iter < 15)
    {
        // Note: The reaction rate is __not__ updated in each iteration!
        auto const qR_updated =
            _react_kin->getReactionRate(p_V, p, T, *solid_state);

        qR = 0.5 * (qR + qR_updated);

        // TODO fix.
        rhoSR[0] = rhoSR_prev_ts + delta_t * qR[0];
    }

    _damped = _damped || reaction_damped;
    return reaction_damped;
}

bool ReactionRateUpdateOnce::isStateCompatible(ReactiveSolidState& state) const
{
    return dynamic_cast<ReactiveSolidStateOneComponent*>(&state) != nullptr;
}

void ReactionRateUpdateOnce::preIteration(const unsigned iter)
{
    _iter = iter;
    // _damped_previous_iteration = _damped;
    _damped = false;
}

void ReactionRateUpdateOnce::postIteration()
{
    if (_damped)
    {
        INFO(
            "ReactionRateUpdateOnce: The reaction rate has been damped for "
            "some elements in the iteration just finished.");
        _last_damped_iter = _iter;
    }
}

std::unique_ptr<ReactionRateData>
ReactionRateUpdateOnce::createReactionRateData() const
{
    return std::unique_ptr<ReactionRateData>(new ReactionRateDataUpdateOnce);
}

std::unique_ptr<ReactionRateUpdateOnce> createReactionRateUpdateOnce(
    BaseLib::ConfigTree const& config)
{
    config.checkConfigParameter("type", "UpdateOnce");

    auto const eq_react_lim = config.getConfigParameter<double>(
        "equilibrium_reaction_vapour_pressure_limit");

    auto const hard_pV_lim =
        config.getConfigParameter<double>("hard_lower_vapour_pressure_limit");

    auto const upper_pV_lim =
        config.getConfigParameter<double>("upper_vapour_pressure_limit");

    auto react_kin =
        createReactionKinetics(config.getConfigSubtree("reaction_kinetics"));

    return std::unique_ptr<ReactionRateUpdateOnce>(new ReactionRateUpdateOnce(
        eq_react_lim, hard_pV_lim, upper_pV_lim, std::move(react_kin)));
}

}  // namespace MaterialLib
