/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ReactionRateUpdateMultipleEnerchem.h"

#include <cassert>
#include <logog/include/logog.hpp>

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "BaseLib/Functional.h"
#include "MaterialLib/PhysicalConstant.h"
#include "MathLib/Nonlinear/Root1D.h"

#include "MaterialLib/ReactionKinetics/ReactionRateDataUpdateMultiple.h"
#include "MaterialLib/ReactionKinetics/ReactiveSolidStateTwoComponents.h"

namespace MaterialLib
{
double ReactionRateUpdateMultipleEnerchem::computeEquilibriumVapourPressure(
    const double p, const double T, ReactiveSolidState const& solid_state) const
{
    double p_V;

    // For low vapour mass fraction compute reaction equilibrium.
    auto const f = [this, p, T, &solid_state](double curr_p_V) -> double {
        auto qR = _react_kin->getReactionRate(curr_p_V, p, T, solid_state);
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

bool ReactionRateUpdateMultipleEnerchem::computeReactionRate(
    double const delta_t, double const p, double const T, double& p_V,
    ReactiveSolidRate& qR, ReactiveSolidState* solid_state,
    ReactionRateData* reaction_rate_data)
{
    const double RATE_FACTOR = 1e-3;

    if (qR.size() != 2)
    {
        OGS_FATAL("error");
    }
    auto const qR_sum = qR.sum();

    auto& conv = solid_state->conversion();
    if (conv.size() != 2)
    {
        OGS_FATAL("error");
    }

    // INFO("qR0 = %e, conv0 = %e.", qR[0], conv[0]);

    assert(dynamic_cast<ReactionRateDataUpdateMultiple*>(reaction_rate_data) !=
           nullptr);
    auto& d = *static_cast<ReactionRateDataUpdateMultiple*>(reaction_rate_data);

    bool reaction_damped = false;

    if (p_V <= 0.0 ||
        (p_V < _equilibrium_reaction_vapour_pressure_limit
         // TODO make configurable
         && std::abs(qR_sum) > 1e-3 * RATE_FACTOR) ||
        p_V > _upper_vapour_pressure_limit)
    {
        // TODO that's new
        // assert(std::abs(qR) <= std::abs(d.reaction_rate_prev_ts[int_pt]));
        d.reaction_rate_max = qR;

        p_V = computeEquilibriumVapourPressure(p, T, *solid_state);

        if (std::abs(qR_sum) > 2e-5 * RATE_FACTOR)
        {
            // Can be interpreted as average of qR at the beginning of the
            // timestep
            // and 0.0 (i.e. equilibrium reached) at the end of the timestep.
            qR /= 2.0;
        }
        else if (qR_sum > -1e-5 * RATE_FACTOR)
        {
            qR.setConstant(-5e-6 * RATE_FACTOR);
        }

        reaction_damped = true;
    }
    else if (_iter > _maximum_reaction_rate_updates)
    {
        return reaction_damped;
    }
    else if (!d.damped_in_current_ts)
    {
        auto const qR_prev_iter = qR;
        auto const qR_curr_iter =
            _react_kin->getReactionRate(p_V, p, T, *solid_state);
        assert(qR_curr_iter.size() == 2);

        // this is an average value over all iterations
        qR = (_iter * qR_prev_iter + qR_curr_iter) / (_iter + 1);
    }
    else
    {
        auto const qR_prev_iter = qR;
        auto const qR_max = d.reaction_rate_max;

        // try to get closer to qR_max (with lessening intensity)
        qR = (_iter * qR_prev_iter + qR_max) / (_iter + 1);
    }

    conv = d.solid_density_prev_ts + delta_t * qR;
    // INFO("--> qR0 = %e, conv0 = %e.", qR[0], conv[0]);

    d.damped_in_current_ts = d.damped_in_current_ts || reaction_damped;
    _damped = _damped || reaction_damped;
    return reaction_damped;
}

bool ReactionRateUpdateMultipleEnerchem::isStateCompatible(
    ReactiveSolidState& state) const
{
    return dynamic_cast<ReactiveSolidStateTwoComponents*>(&state) != nullptr;
}

void ReactionRateUpdateMultipleEnerchem::preIteration(const unsigned iter)
{
    INFO("ReactionRateUpdateMultipleEnerchem pre iter");
    _iter = iter;
    // _damped_previous_iteration = _damped;
    _damped = false;
}

void ReactionRateUpdateMultipleEnerchem::postIteration()
{
    if (_damped)
    {
        INFO(
            "ReactionRateUpdateMultipleEnerchem: The reaction rate has been "
            "damped for "
            "some elements in the iteration just finished.");
        _last_damped_iter = _iter;
    }
}

std::vector<NumLib::NamedFunction>
ReactionRateUpdateMultipleEnerchem::getNamedFunctions() const
{
    /*
    auto fcts = _react_kin->getNamedFunctions();
    auto cb = [&](const double p_V, const double p, const double T,
                  const double rho_SR) -> double {
        ReactiveSolidStateTwoComponents s{{rho_SR}};
        auto qR = _react_kin->getReactionRate(p_V, p, T, s);
        assert(qR.size() == 1);
        return qR[0];
    };

    fcts.emplace_back(NumLib::NamedFunction{
        "raw_reaction_rate",
        {"vapour_partial_pressure", "pressure", "temperature", "solid_density"},
        BaseLib::easyBind(std::move(cb))});
    */

    return _react_kin->getNamedFunctions();
}

std::unique_ptr<ReactionRateData>
ReactionRateUpdateMultipleEnerchem::createReactionRateData() const
{
    return std::unique_ptr<ReactionRateData>(
        new ReactionRateDataUpdateMultiple);
}

std::unique_ptr<ReactionRateUpdateMultipleEnerchem>
createReactionRateUpdateMultipleEnerchem(BaseLib::ConfigTree const& config)
{
    config.checkConfigParameter("type", "UpdateMultipleEnerchem");

    auto const eq_react_lim = config.getConfigParameter<double>(
        "equilibrium_reaction_vapour_pressure_limit");

    auto const max_updates =
        config.getConfigParameter<unsigned>("maximum_reaction_rate_updates");

    auto const upper_pV_lim =
        config.getConfigParameter<double>("upper_vapour_pressure_limit");

    auto react_kin =
        createReactionKinetics(config.getConfigSubtree("reaction_kinetics"));

    return std::unique_ptr<ReactionRateUpdateMultipleEnerchem>(
        new ReactionRateUpdateMultipleEnerchem(
            eq_react_lim, upper_pV_lim, max_updates, std::move(react_kin)));
}

}  // namespace MaterialLib
