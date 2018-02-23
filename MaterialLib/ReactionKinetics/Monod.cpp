/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Monod.h"
#include "BaseLib/ConfigTree.h"
#include "BaseLib/Functional.h"

namespace MaterialLib
{
std::vector<NumLib::NamedFunction> MonodKinetics::getNamedFunctions() const
{
    auto fcts = _raw_kinetics->getNamedFunctions();
    fcts.emplace_back(NumLib::NamedFunction{
        "monod_factor",
        {"pressure"},
        BaseLib::easyBind(&MonodKinetics::getMonodTerm, this)});

    return fcts;
}

std::unique_ptr<ReactionKinetics> createMonodKinetics(
    BaseLib::ConfigTree const& config)
{
    config.checkConfigParameter("type", "Monod");

    auto const p_half = config.getConfigParameter<double>("p_half");
    auto raw_kinetics =
        createReactionKinetics(config.getConfigSubtree("raw_kinetics"));

    return std::unique_ptr<ReactionKinetics>(
        new MonodKinetics(p_half, std::move(raw_kinetics)));
}

}  // namespace MaterialLib
