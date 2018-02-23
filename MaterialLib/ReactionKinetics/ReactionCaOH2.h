/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ReactionKinetics.h"

namespace MaterialLib
{
class ReactionCaOH2 final : public ReactionKinetics
{
public:
    ReactiveSolidRate getReactionRate(
        const double p_V, const double /*p*/, const double T,
        ReactiveSolidState const& solid_state) override;

    HeatOfReactionData getHeatOfReaction(
        const double p_V, const double T,
        const ReactiveSolidState* const state) const override;

    bool isStateCompatible(ReactiveSolidState& state) const override
    {
        return state.conversion().size() == 1;
    }

private:
    //! reaction enthalpy in J/mol; negative for exothermic composition reaction
    static const double _reaction_enthalpy;
    static const double _reaction_entropy;  //!< reaction entropy in J/mol/K
    static const double _M_react;           //!< reactive component molar mass

    static const double _tol_l;
    static const double _tol_u;
    static const double _tol_rho;

public:
    static const double rho_low;  //! lower density limit
    static const double rho_up;   //! upper density limit
};

std::unique_ptr<ReactionCaOH2> createReactionCaOH2(
    BaseLib::ConfigTree const& config);

}  // namespace MaterialLib
