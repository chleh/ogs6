/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LDFTwoComponentTest3_LoadingPV.h"

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/ReactionKinetics/Enerchem/SorptionEquilibriumCaCl2CaX_2_SeparateLoadingPV.h"
#include "MaterialLib/ReactionKinetics/ReactiveSolidStateTwoComponents.h"

namespace MaterialLib
{
LDFTwoComponentTest3_LoadingPV::LDFTwoComponentTest3_LoadingPV(
    std::unique_ptr<ReactionEquilibrium2>&& equil,
    std::vector<double> const& k_LDF_)
    : _equil(std::move(equil)), _k_LDF(MathLib::toVector(k_LDF_))
{
    if (_k_LDF.size() != 2)
        OGS_FATAL("k_LDF has wrong number of components");

    if (!dynamic_cast<SorptionEquilibriumCaCl2CaX_2_SeparateLoadingPV*>(
            _equil.get()))
        OGS_FATAL("wrong reaction equilibrium");
}

ReactiveSolidRate LDFTwoComponentTest3_LoadingPV::getReactionRate(
    const double p_V, const double /*p*/, const double T,
    const ReactiveSolidState& solid_state)
{
    ReactiveSolidStateTwoComponents eq{{0.0, 0.0}};
    _equil->getEquilibrium(p_V, T, solid_state, eq);
    return _k_LDF.cwiseProduct(eq.conversion() - solid_state.conversion());
}

HeatOfReactionData LDFTwoComponentTest3_LoadingPV::getHeatOfReaction(
    const double p_V, const double T,
    const ReactiveSolidState* const state) const
{
    return _equil->getHeatOfReaction(p_V, T, state);
}

std::unique_ptr<LDFTwoComponentTest3_LoadingPV>
createLDFTwoComponentTest3_LoadingPV(BaseLib::ConfigTree const& config)
{
    config.checkConfigParameter("type", "LDFTwoComponentTest3_LoadingPV");

    auto equil =
        createReactionEquilibrium2(config.getConfigSubtree("equilibrium"));

    auto const k_LDF = config.getConfigParameter<std::vector<double>>("k_ldf");

    return std::unique_ptr<LDFTwoComponentTest3_LoadingPV>{
        new LDFTwoComponentTest3_LoadingPV{std::move(equil), k_LDF}};
}

}  // namespace MaterialLib
