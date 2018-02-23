/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MaterialLib/ReactionKinetics/ReactionKinetics.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"

#include "ReactionEquilibrium2.h"

namespace MaterialLib
{
class LDFTwoComponentTest4_LoadingLoadingDamped : public ReactionKinetics
{
public:
    explicit LDFTwoComponentTest4_LoadingLoadingDamped(
        std::unique_ptr<ReactionEquilibrium2>&& equil,
        std::vector<double> const& k_LDF_,
        const double zeolite_kinetics_damping_factor);
    ReactiveSolidRate getReactionRate(
        const double p_V, const double /*p*/, const double T,
        ReactiveSolidState const& solid_state) override;

    HeatOfReactionData getHeatOfReaction(
        const double p_V, const double T,
        ReactiveSolidState const* const state) const override;

    std::vector<NumLib::NamedFunction> getNamedFunctions() const override
    {
        return _equil->getNamedFunctions();
    }

    bool isStateCompatible(ReactiveSolidState& state) const override
    {
        return state.conversion().size() == 2;
    }

private:
    std::unique_ptr<ReactionEquilibrium2> _equil;
    Eigen::VectorXd const _k_LDF;
    const double _zeolite_kinetics_damping_factor;
};

std::unique_ptr<LDFTwoComponentTest4_LoadingLoadingDamped>
createLDFTwoComponentTest4_LoadingLoadingDamped(
    BaseLib::ConfigTree const& config);

}  // namespace MaterialLib
