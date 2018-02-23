/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LinearDrivingForceGorbachCoefficient.h"
#include "BaseLib/ConfigTree.h"
#include "ReactionEquilibrium.h"

namespace MaterialLib
{
std::unique_ptr<ReactionKinetics> createLinearDrivingForceGorbachCoefficient(
    BaseLib::ConfigTree const& config)
{
    config.checkConfigParameter("type", "LinearDrivingForceGorbachCoefficient");

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

    auto equil =
        createReactionEquilibrium(config.getConfigSubtree("equilibrium"));

    return std::unique_ptr<ReactionKinetics>(
        new LinearDrivingForceGorbachCoefficient(
            pellet_radius, pellet_porosity, pellet_tortuosity,
            pellet_density_dry, molar_mass_reactive_component,
            std::move(diffusion_coefficient), std::move(equil)));
}

}  // namespace MaterialLib
