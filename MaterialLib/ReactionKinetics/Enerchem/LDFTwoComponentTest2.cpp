/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LDFTwoComponentTest2.h"

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/ReactionKinetics/Enerchem/SorptionEquilibriumCaCl2CaX_1_SeparateLoadingLoading.h"
#include "MaterialLib/ReactionKinetics/ReactiveSolidStateTwoComponents.h"

namespace MaterialLib
{
LDFTwoComponentTest2::LDFTwoComponentTest2(
    std::unique_ptr<ReactionEquilibrium2>&& equil,
    std::vector<double> const& k_LDF_)
    : _equil(std::move(equil)), _k_LDF(MathLib::toVector(k_LDF_))
{
    if (_k_LDF.size() != 2)
        OGS_FATAL("k_LDF has wrong number of components");

    if (!dynamic_cast<SorptionEquilibriumCaCl2CaX_1_SeparateLoadingLoading*>(
            _equil.get()))
        OGS_FATAL("wrong reaction equilibrium");
}

ReactiveSolidRate LDFTwoComponentTest2::getReactionRate(
    const double p_V, const double /*p*/, const double T,
    const ReactiveSolidState& solid_state)
{
    ReactiveSolidStateTwoComponents eq{{0.0, 0.0}};
    _equil->getEquilibrium(p_V, T, solid_state, eq);
    return _k_LDF.cwiseProduct(eq.conversion() - solid_state.conversion());
}

HeatOfReactionData LDFTwoComponentTest2::getHeatOfReaction(
    const double p_V, const double T,
    const ReactiveSolidState* const state) const
{
    return _equil->getHeatOfReaction(p_V, T, state);
}

std::unique_ptr<LDFTwoComponentTest2> createLDFTwoComponentTest2(
    BaseLib::ConfigTree const& config)
{
    config.checkConfigParameter("type", "LDFTwoComponentTest2");

    auto equil =
        createReactionEquilibrium2(config.getConfigSubtree("equilibrium"));

    auto const k_LDF = config.getConfigParameter<std::vector<double>>("k_ldf");

    return std::unique_ptr<LDFTwoComponentTest2>{
        new LDFTwoComponentTest2{std::move(equil), k_LDF}};
}

}  // namespace MaterialLib
