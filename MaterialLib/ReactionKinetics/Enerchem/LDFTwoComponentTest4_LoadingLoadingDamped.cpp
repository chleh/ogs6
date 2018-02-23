/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LDFTwoComponentTest4_LoadingLoadingDamped.h"

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/ReactionKinetics/ReactiveSolidStateTwoComponents.h"

namespace MaterialLib
{
LDFTwoComponentTest4_LoadingLoadingDamped::
    LDFTwoComponentTest4_LoadingLoadingDamped(
        std::unique_ptr<ReactionEquilibrium2>&& equil,
        std::vector<double> const& k_LDF_,
        const double zeolite_kinetics_damping_factor)
    : _equil(std::move(equil)),
      _k_LDF(MathLib::toVector(k_LDF_)),
      _zeolite_kinetics_damping_factor(zeolite_kinetics_damping_factor)
{
    if (_k_LDF.size() != 2)
        OGS_FATAL("k_LDF has wrong number of components");

    if (!dynamic_cast<ReactionEquilibrium2_LoadingLoading*>(_equil.get()))
        OGS_FATAL("wrong reaction equilibrium");
}

ReactiveSolidRate LDFTwoComponentTest4_LoadingLoadingDamped::getReactionRate(
    const double p_V, const double /*p*/, const double T,
    const ReactiveSolidState& solid_state)
{
    ReactiveSolidStateTwoComponents eq{{0.0, 0.0}};
    _equil->getEquilibrium(p_V, T, solid_state, eq);

    auto const C_zeo = solid_state.conversion()[0];

    auto const loading_dependent_damping =
        std::max(1e-3, 1.0 - _zeolite_kinetics_damping_factor * C_zeo);

    auto rate =
        (_k_LDF.cwiseProduct(eq.conversion() - solid_state.conversion()))
            .eval();
    rate[0] *= loading_dependent_damping;
    return rate;
}

HeatOfReactionData LDFTwoComponentTest4_LoadingLoadingDamped::getHeatOfReaction(
    const double p_V, const double T,
    const ReactiveSolidState* const state) const
{
    return _equil->getHeatOfReaction(p_V, T, state);
}

std::unique_ptr<LDFTwoComponentTest4_LoadingLoadingDamped>
createLDFTwoComponentTest4_LoadingLoadingDamped(
    BaseLib::ConfigTree const& config)
{
    config.checkConfigParameter("type",
                                "LDFTwoComponentTest4_LoadingLoadingDamped");

    auto equil =
        createReactionEquilibrium2(config.getConfigSubtree("equilibrium"));

    auto const k_LDF = config.getConfigParameter<std::vector<double>>("k_ldf");

    auto const salt_loading = config.getConfigParameter<double>("salt_loading");

    // Cf.
    // /home/lehmannc/docs/2016-enerchem/data-repo/simulation/05-adsorption-homog-C0-LDF-model-test1-fixed-t-end/opt-params/fit-by-phase.png

    // m:   0.16319249 +/- 0.007794 (4.78%)
    auto const damping_factor = 0.16319249 * salt_loading;

    return std::unique_ptr<LDFTwoComponentTest4_LoadingLoadingDamped>{
        new LDFTwoComponentTest4_LoadingLoadingDamped{std::move(equil), k_LDF,
                                                      damping_factor}};
}

}  // namespace MaterialLib
