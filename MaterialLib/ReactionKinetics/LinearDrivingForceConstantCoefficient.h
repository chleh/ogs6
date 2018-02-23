/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#pragma once

#include "ReactionEquilibrium.h"
#include "ReactionKineticsByEquilibrium.h"

namespace MaterialLib
{
class LinearDrivingForceConstantCoefficient final
    : public ReactionKineticsByEquilibrium
{
public:
    LinearDrivingForceConstantCoefficient(
        const double k_rate, std::unique_ptr<ReactionEquilibrium>&& equil)
        : _k_rate(k_rate), _equil(std::move(equil))
    {
    }

    ReactiveSolidRate getReactionRate(
        const double p_V, const double /*p*/, const double T,
        ReactiveSolidState const& solid_state) override
    {
        assert(solid_state.conversion().size() == 1);
        return Eigen::VectorXd::Constant(
            1, _k_rate * (_equil->getEquilibriumDensity(p_V, T) -
                          solid_state.conversion()[0]));
    }

    HeatOfReactionData getHeatOfReaction(
        const double p_V, const double T,
        const ReactiveSolidState* const state) const override
    {
        return _equil->getHeatOfReaction(p_V, T, state);
    }

    std::vector<NumLib::NamedFunction> getNamedFunctions() const override
    {
        return _equil->getNamedFunctions();
    }

    bool isStateCompatible(ReactiveSolidState& state) const override
    {
        return state.conversion().size() == 1;
    }

private:
    const double _k_rate;
    std::unique_ptr<ReactionEquilibrium> _equil;
};

std::unique_ptr<ReactionKinetics> createLinearDrivingForceConstantCoefficient(
    BaseLib::ConfigTree const& config);

}  // namespace MaterialLib
