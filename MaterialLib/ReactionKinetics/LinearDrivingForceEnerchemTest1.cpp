/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LinearDrivingForceEnerchemTest1.h"
#include "BaseLib/ConfigTree.h"
#include "ReactionEquilibrium.h"

namespace MaterialLib
{
std::unique_ptr<ReactionKinetics> createLinearDrivingForceEnerchemTest1(
    BaseLib::ConfigTree const& config)
{
    config.checkConfigParameter("type", "LinearDrivingForceEnerchemTest1");

    auto const pellet_radius =
        config.getConfigParameter<double>("pellet_radius");
    auto const pellet_porosity =
        config.getConfigParameter<double>("pellet_porosity");
    auto const pellet_tortuosity =
        config.getConfigParameter<double>("pellet_tortuosity");
    auto const pellet_density_dry =
        config.getConfigParameter<double>("pellet_density_dry");
    auto const molar_mass_reactive_component =
        config.getConfigParameter<double>("molar_mass_reactive_component");
    auto diffusion_coefficient = ProcessLib::TES::createDiffusionCoefficient(
        config.getConfigSubtree("diffusion_coefficient"));

    auto const loading_dependent_damping_factor =
        config.getConfigParameterOptional<double>(
            "loading_dependent_damping_factor");

    double df;
    if (loading_dependent_damping_factor)
        df = *loading_dependent_damping_factor;
    else
    {
        auto const salt_loading =
            config.getConfigParameter<double>("salt_loading");

        // Cf.
        // /home/lehmannc/docs/2016-enerchem/data-repo/simulation/05-adsorption-homog-C0-LDF-model-test1-fixed-t-end/opt-params/fit-by-phase.png

        // m:   0.16319249 +/- 0.007794 (4.78%)
        df = 0.16319249 * salt_loading;
    }

    auto equil =
        createReactionEquilibrium(config.getConfigSubtree("equilibrium"));

    return std::unique_ptr<ReactionKinetics>(
        new LinearDrivingForceEnerchemTest1(
            pellet_radius, pellet_porosity, pellet_tortuosity,
            pellet_density_dry, df, molar_mass_reactive_component,
            std::move(diffusion_coefficient), std::move(equil)));
}

}  // namespace MaterialLib
