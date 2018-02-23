/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ReactionEquilibrium.h"
#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"

#include "DubininPolanyi.h"
#include "DubininPolanyiOrTPD.h"
#include "Enerchem/SorptionEquilibriumCaCl2CaX_5_DubDub.h"
#include "Enerchem/SorptionEquilibriumCaCl2CaX_6_DubDub.h"
#include "SorptionEquilibriumCaCl2CaX_1_Combined.h"
#include "SorptionEquilibriumCaCl2CaX_1_Combined_ZeoHeat.h"

namespace MaterialLib
{
std::unique_ptr<ReactionEquilibrium> createReactionEquilibrium(
    BaseLib::ConfigTree const& config)
{
    auto const type = config.peekConfigParameter<std::string>("type");

    if (type == "DubininPolanyi")
        return createDubininPolanyi(config);
    else if (type == "DubininPolanyiOrTPD")
        return createDubininPolanyiOrTPD(config);
    else if (type == "SorptionEquilibriumCaCl2CaX_1_Combined")
        return createSorptionEquilibriumCaCl2CaX_1_Combined(config);
    else if (type == "SorptionEquilibriumCaCl2CaX_1_Combined_ZeoHeat")
        return createSorptionEquilibriumCaCl2CaX_1_Combined_ZeoHeat(config);
    else if (type == "SorptionEquilibriumCaCl2CaX_5_DubDub")
        return createSorptionEquilibriumCaCl2CaX_5_DubDub(config);
    else if (type == "SorptionEquilibriumCaCl2CaX_6_DubDub")
        return createSorptionEquilibriumCaCl2CaX_6_DubDub(config);

    OGS_FATAL("There is no reaction equilibrium of type `%s'.", type.c_str());
}

}  // namespace MaterialLib
