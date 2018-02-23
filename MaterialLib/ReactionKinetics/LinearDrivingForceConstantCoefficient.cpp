/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LinearDrivingForceConstantCoefficient.h"
#include "BaseLib/ConfigTree.h"
#include "ReactionEquilibrium.h"

namespace MaterialLib
{
std::unique_ptr<ReactionKinetics> createLinearDrivingForceConstantCoefficient(
    BaseLib::ConfigTree const& config)
{
    config.checkConfigParameter("type",
                                "LinearDrivingForceConstantCoefficient");

    auto const k_rate = config.getConfigParameter<double>("k_rate");
    auto equil =
        createReactionEquilibrium(config.getConfigSubtree("equilibrium"));

    return std::unique_ptr<ReactionKinetics>(
        new LinearDrivingForceConstantCoefficient(k_rate, std::move(equil)));
}

}  // namespace MaterialLib
