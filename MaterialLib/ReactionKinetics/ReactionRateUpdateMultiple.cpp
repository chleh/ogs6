/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ReactionRateUpdateMultiple.h"

#include <cassert>
#include <logog/include/logog.hpp>

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "BaseLib/Functional.h"
#include "MaterialLib/PhysicalConstant.h"
#include "MathLib/Nonlinear/Root1D.h"

#include "ReactionRateDataUpdateMultiple.h"
#include "ReactiveSolidStateOneComponent.h"

namespace MaterialLib
{
double ReactionRateUpdateMultiple::computeEquilibriumVapourPressure(
    const double p, const double T, const double rhoSR) const
{
    double p_V;

    // For low vapour mass fraction compute reaction equilibrium.
    auto const f = [this, p, T, rhoSR](double curr_p_V) -> double {
        ReactiveSolidStateOneComponent s{{rhoSR}};
        auto qR = _react_kin->getReactionRate(curr_p_V, p, T, s);
        assert(qR.size() == 1);
        return qR[0];
    };
    // Find root of reaction rate (== reaction equilibrium).
    auto rf = MathLib::Nonlinear::makeRegulaFalsi<MathLib::Nonlinear::Pegasus>(
        f, 1e-6, _upper_vapour_pressure_limit);

    // Note: Modifying p_V deteriorates mass conservation!
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
            p_V = 1e-6;
        /*
            OGS_FATAL(
                "Regula Falsi failed. The signs of function values at both "
                "interval ends are the same."
                " The smaller function absolute value is %g at %g",
                p_q.second, p_q.first);
                */
    }

    return p_V;
}

bool ReactionRateUpdateMultiple::computeReactionRate(
    double const delta_t, double const p, double const T, double& p_V,
    ReactiveSolidRate& qR, ReactiveSolidState* solid_state,
    ReactionRateData* reaction_rate_data)
{
    assert(qR.size() == 1);
    auto& qR0 = qR[0];

    auto& rhoSR = solid_state->conversion();
    assert(rhoSR.size() == 1);

    assert(dynamic_cast<ReactionRateDataUpdateMultiple*>(reaction_rate_data) !=
           nullptr);
    auto& d = *static_cast<ReactionRateDataUpdateMultiple*>(reaction_rate_data);
    auto const rhoSR_prev_ts = d.solid_density_prev_ts[0];
    auto const damped_in_current_ts = d.damped_in_current_ts;

    bool reaction_damped = false;

    if (p_V <= 0.0 ||
        (p_V < _equilibrium_reaction_vapour_pressure_limit
         // TODO make configurable
         && std::abs(qR0) > 1e-3) ||
        p_V > _upper_vapour_pressure_limit)
    {
        // TODO that's new
        // assert(std::abs(qR) <= std::abs(d.reaction_rate_prev_ts[int_pt]));
        d.reaction_rate_max.setConstant(qR0);

        p_V = computeEquilibriumVapourPressure(p, T, rhoSR[0]);

        if (std::abs(qR0) > 2e-5)
        {
            // Can be interpreted as average of qR at the beginning of the
            // timestep
            // and 0.0 (i.e. equilibrium reached) at the end of the timestep.
            qR0 /= 2.0;
        }
        else if (qR0 > -1e-5)
        {
            qR0 = -1e-5;
        }

        reaction_damped = true;
    }
    else if (_iter > _maximum_reaction_rate_updates)
    {
        return reaction_damped;
    }
    else if (!damped_in_current_ts)
    {
        auto const qR_prev_iter = qR0;
        auto const qR_curr_iter =
            _react_kin->getReactionRate(p_V, p, T, *solid_state);
        assert(qR_curr_iter.size() == 1);

        // this is an average value over all iterations
        qR0 = (_iter * qR_prev_iter + qR_curr_iter[0]) / (_iter + 1);
    }
    else
    {
        auto const qR_prev_iter = qR0;
        auto const qR_max = d.reaction_rate_max[0];

        // try to get closer to qR_max (with lessening intensity)
        qR0 = (_iter * qR_prev_iter + qR_max) / (_iter + 1);
    }

    rhoSR[0] = rhoSR_prev_ts + delta_t * qR0;

    d.damped_in_current_ts = damped_in_current_ts || reaction_damped;
    _damped = _damped || reaction_damped;
    return reaction_damped;
}

bool ReactionRateUpdateMultiple::isStateCompatible(
    ReactiveSolidState& state) const
{
    return dynamic_cast<ReactiveSolidStateOneComponent*>(&state) != nullptr;
}

void ReactionRateUpdateMultiple::preIteration(const unsigned iter)
{
    _iter = iter;
    // _damped_previous_iteration = _damped;
    _damped = false;
}

void ReactionRateUpdateMultiple::postIteration()
{
    if (_damped)
    {
        INFO(
            "ReactionRateUpdateMultiple: The reaction rate has been damped for "
            "some elements in the iteration just finished.");
        _last_damped_iter = _iter;
    }
}

std::vector<NumLib::NamedFunction>
ReactionRateUpdateMultiple::getNamedFunctions() const
{
    auto fcts = _react_kin->getNamedFunctions();
    auto cb = [&](const double p_V, const double p, const double T,
                  const double rho_SR) -> double {
        ReactiveSolidStateOneComponent s{{rho_SR}};
        auto qR = _react_kin->getReactionRate(p_V, p, T, s);
        assert(qR.size() == 1);
        return qR[0];
    };

    fcts.emplace_back(NumLib::NamedFunction{
        "raw_reaction_rate",
        {"vapour_partial_pressure", "pressure", "temperature", "solid_density"},
        BaseLib::easyBind(std::move(cb))});

    return fcts;
}

std::unique_ptr<ReactionRateData>
ReactionRateUpdateMultiple::createReactionRateData() const
{
    return std::unique_ptr<ReactionRateData>(
        new ReactionRateDataUpdateMultiple);
}

std::unique_ptr<ReactionRateUpdateMultiple> createReactionRateUpdateMultiple(
    BaseLib::ConfigTree const& config)
{
    config.checkConfigParameter("type", "UpdateMultiple");

    auto const eq_react_lim = config.getConfigParameter<double>(
        "equilibrium_reaction_vapour_pressure_limit");

    auto const max_updates =
        config.getConfigParameter<unsigned>("maximum_reaction_rate_updates");

    auto const upper_pV_lim =
        config.getConfigParameter<double>("upper_vapour_pressure_limit");

    auto react_kin =
        createReactionKinetics(config.getConfigSubtree("reaction_kinetics"));

    return std::unique_ptr<ReactionRateUpdateMultiple>(
        new ReactionRateUpdateMultiple(eq_react_lim, upper_pV_lim, max_updates,
                                       std::move(react_kin)));
}

}  // namespace MaterialLib
