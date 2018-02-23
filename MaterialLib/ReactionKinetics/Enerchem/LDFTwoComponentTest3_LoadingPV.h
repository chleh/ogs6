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
class LDFTwoComponentTest3_LoadingPV : public ReactionKinetics
{
public:
    explicit LDFTwoComponentTest3_LoadingPV(
        std::unique_ptr<ReactionEquilibrium2>&& equil,
        std::vector<double> const& k_LDF_);
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
};

std::unique_ptr<LDFTwoComponentTest3_LoadingPV>
createLDFTwoComponentTest3_LoadingPV(BaseLib::ConfigTree const& config);

}  // namespace MaterialLib
