/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#pragma once

#include "ReactionEquilibrium2.h"

namespace MaterialLib
{
/* Dasselbe wie SorptionEquilibriumCaCl2CaX_1_SeparateLoadingLoading,
 * nur dass die Auswaschung des Salzes aus den Zeolithporen ignoriert wird:
 * Es bildet sich eine Salzlösung in den Makroporen, aber der Zeolith bleibt
 * blockiert. Mglw. kann damit der Knick in der Durchbruchskurve verringert
 * werden.
 */
class SorptionEquilibriumCaCl2CaX_3_SeparateLoadingLoading
    : public ReactionEquilibrium2_LoadingLoading
{
public:
    SorptionEquilibriumCaCl2CaX_3_SeparateLoadingLoading(
        const double salt_loading)
        : _salt_loading(salt_loading)
    {
    }

    void getEquilibrium(const double p_Ads,
                        const double T_Ads,
                        ReactiveSolidState const& solid_state,
                        ReactiveSolidState& equilibrium) const override;

    HeatOfReactionData getHeatOfReaction(
        const double p_Ads, const double T_Ads,
        ReactiveSolidState const* const state) const override;

    std::vector<NumLib::NamedFunction> getNamedFunctions() const override;

protected:
    const double _salt_loading;  //!< in mol/mol
};

std::unique_ptr<SorptionEquilibriumCaCl2CaX_3_SeparateLoadingLoading>
createSorptionEquilibriumCaCl2CaX_3_SeparateLoadingLoading(
    BaseLib::ConfigTree const& config);

}  // namespace MaterialLib
