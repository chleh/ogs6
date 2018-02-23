/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ReactiveSolidModelDensityBased.h"
#include "ReactiveSolidModelLoadingBased.h"

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"

namespace MaterialLib
{
std::unique_ptr<ReactiveSolidModel> createReactiveSolidModel(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters)
{
    auto const type = config.peekConfigParameter<std::string>("type");

    if (type == "LoadingBased")
        return createReactiveSolidModelLoadingBased(config, parameters);
    else if (type == "DensityBased")
        return createReactiveSolidModelDensityBased(config, parameters);

    OGS_FATAL("Unknown reactive solid type `%s'", type.c_str());
}

}  // namespace MaterialLib
