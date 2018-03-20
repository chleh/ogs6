#pragma once

#include <unordered_map>

#include "BaseLib/ConfigTree.h"

#include "TCHSNoStokesMaterial.h"

namespace ProcessLib
{
namespace TCHSNoStokes
{
namespace Material
{
TCHSNoStokesMaterial createTCHSNoStokesMaterial(
    BaseLib::ConfigTree const& config,
    const std::vector<std::unique_ptr<ParameterBase>>& parameters);

std::unordered_map<int, TCHSNoStokesMaterial> createTCHSNoStokesMaterials(
    BaseLib::ConfigTree const& config,
    const std::vector<std::unique_ptr<ParameterBase>>& parameters);

std::unique_ptr<Permeability> createPermeability(
    BaseLib::ConfigTree const& config);

}  // namespace Material

}  // namespace TCHSNoStokes

}  // namespace ProcessLib
