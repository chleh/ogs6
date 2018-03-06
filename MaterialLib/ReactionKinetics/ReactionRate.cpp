/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ReactionRate.h"

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"

#include "Enerchem/ReactionRateUpdateMultipleEnerchem.h"
#include "ReactionRateConstant.h"
#include "ReactionRateInert.h"
#include "ReactionRateRaw.h"
#include "ReactionRateRawWithODESolver.h"
#include "ReactionRateSimple.h"
#include "ReactionRateUpdateMultiple.h"
#include "ReactionRateUpdateOnce.h"

namespace MaterialLib
{
std::unique_ptr<ReactionRate> createReactionRate(
    const BaseLib::ConfigTree& config)
{
    auto const type = config.peekConfigParameter<std::string>("type");

    if (type == "UpdateOnce")
    {
        return createReactionRateUpdateOnce(config);
    }
    else if (type == "UpdateMultiple")
    {
        return createReactionRateUpdateMultiple(config);
    }
    else if (type == "Raw")
    {
        return createReactionRateRaw(config);
    }
    else if (type == "Simple")
    {
        return createReactionRateSimple(config);
    }
    else if (type == "Inert")
    {
        config.ignoreConfigParameter("type");
        return std::unique_ptr<ReactionRate>(new ReactionRateInert);
    }
    else if (type == "Constant")
    {
        return createReactionRateConstant(config);
    }
    else if (type == "ReactionRateRawWithODESolver")
    {
        return createReactionRateRawWithODESolver(config);
    }
    else if (type == "UpdateMultipleEnerchem")
    {
        return createReactionRateUpdateMultipleEnerchem(config);
    }

    OGS_FATAL("There is no reaction rate class of type `%s'.", type.c_str());
}

}  // namespace MaterialLib
