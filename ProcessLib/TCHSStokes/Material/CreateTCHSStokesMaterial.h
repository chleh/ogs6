#pragma once

#include <unordered_map>

#include "BaseLib/ConfigTree.h"

#include "TCHSStokesMaterial.h"

namespace ProcessLib
{
namespace TCHSStokes
{
namespace Material
{
TCHSStokesMaterial createTCHSStokesMaterial(BaseLib::ConfigTree const& config);

std::unordered_map<int, TCHSStokesMaterial> createTCHSStokesMaterials(
    BaseLib::ConfigTree const& config);

std::unique_ptr<FluidDensity> createFluidDensity(
    BaseLib::ConfigTree const& config);

std::unique_ptr<FluidViscosity> createFluidViscosity(
    BaseLib::ConfigTree const& config);

std::unique_ptr<FluidHeatCapacity> createFluidHeatCapacity(
    BaseLib::ConfigTree const& config);

std::unique_ptr<HeatConductivity> createHeatConductivity(
    BaseLib::ConfigTree const& config);

std::unique_ptr<FluidMomentumProductionCoefficient>
createFluidMomentumProductionCoefficient(BaseLib::ConfigTree const& config);

std::unique_ptr<SolidHeatCapacity> createSolidHeatCapacity(
    BaseLib::ConfigTree const& config);

std::unique_ptr<Porosity> createPorosity(BaseLib::ConfigTree const& config);

}  // namespace Material

}  // namespace TCHSStokes

}  // namespace ProcessLib
