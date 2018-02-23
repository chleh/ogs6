/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATERIALLIB_MONODKINETICS_H
#define MATERIALLIB_MONODKINETICS_H

#include "ReactionKinetics.h"

namespace MaterialLib
{
class MonodKinetics final : public ReactionKinetics
{
public:
    MonodKinetics(const double p_half,
                  std::unique_ptr<ReactionKinetics>&& raw_kinetics)
        : _p_half(p_half), _raw_kinetics(std::move(raw_kinetics))
    {
    }

    ReactiveSolidRate getReactionRate(
        const double p_V, const double p, const double T,
        ReactiveSolidState const& solid_state) override
    {
        const double monod_term = p_V / (p_V + _p_half);
        return _raw_kinetics->getReactionRate(p_V, p, T, solid_state) *
               monod_term;
    }

    HeatOfReactionData getHeatOfReaction(
        const double p_V, const double T,
        const ReactiveSolidState* const state) const
    {
        return _raw_kinetics->getHeatOfReaction(p_V, T, state);
    }

    std::vector<NumLib::NamedFunction> getNamedFunctions() const override;

    bool isStateCompatible(ReactiveSolidState& state) const override
    {
        return _raw_kinetics->isStateCompatible(state);
    }

private:
    double getMonodTerm(const double p_V) const
    {
        return p_V / (p_V + _p_half);
    }

    const double _p_half;
    std::unique_ptr<ReactionKinetics> _raw_kinetics;
};

std::unique_ptr<ReactionKinetics> createMonodKinetics(
    BaseLib::ConfigTree const& config);

}  // namespace MaterialLib

#endif  // MATERIALLIB_MONODKINETICS_H
