/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateMFront.h"

#include <dlfcn.h>

namespace MaterialLib
{
namespace Solids
{
namespace MFront
{
template <int DisplacementDim>
std::unique_ptr<MFront<DisplacementDim>> createMFront(
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config)
{
    config.checkConfigParameter("type", "MFront");

    // TODO pass prj file dir.
    auto const lib_path = config.getConfigParameter<std::string>("library");
    auto const model_name = config.getConfigParameter<std::string>("model");

    dlerror();
    auto* lib = dlopen(lib_path.c_str(), RTLD_NOW);
    auto* err = dlerror();
    if (err != nullptr)
    {
        OGS_FATAL("Could not load shared library `%s'. The error is: `%s'.",
                  lib_path.c_str(), err);
    }

    return nullptr;
}

template std::unique_ptr<MFront<2>> createMFront<2>(
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config);
template std::unique_ptr<MFront<3>> createMFront<3>(
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config);

}  // namespace MFront
}  // namespace Solids
}  // namespace MaterialLib
