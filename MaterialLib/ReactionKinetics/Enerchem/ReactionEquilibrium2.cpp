/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ReactionEquilibrium2.h"

#include "BaseLib/ConfigTree.h"
#include "SorptionEquilibriumCaCl2CaX_1_SeparateLoadingLoading.h"
#include "SorptionEquilibriumCaCl2CaX_2_SeparateLoadingPV.h"
#include "SorptionEquilibriumCaCl2CaX_3_SeparateLoadingLoading.h"
#include "SorptionEquilibriumCaCl2CaX_4_SeparateLoadingLoading.h"

namespace MaterialLib
{
std::unique_ptr<ReactionEquilibrium2> createReactionEquilibrium2(
    BaseLib::ConfigTree const& config)
{
    auto const type = config.peekConfigParameter<std::string>("type");

    if (type == "CaCl2CaX_1_SeparateLoadingLoading")
        return createSorptionEquilibriumCaCl2CaX_1_SeparateLoadingLoading(
            config);
    else if (type == "CaCl2CaX_1_SeparateLoadingPV")
        return createSorptionEquilibriumCaCl2CaX_2_SeparateLoadingPV(config);
    else if (type == "CaCl2CaX_3_SeparateLoadingLoading")
        return createSorptionEquilibriumCaCl2CaX_3_SeparateLoadingLoading(
            config);
    else if (type == "CaCl2CaX_4_SeparateLoadingLoading")
        return createSorptionEquilibriumCaCl2CaX_4_SeparateLoadingLoading(
            config);

    OGS_FATAL("Unknown reaction equilibrium: `%s'.", type.c_str());
}

}  // namespace MaterialLib
