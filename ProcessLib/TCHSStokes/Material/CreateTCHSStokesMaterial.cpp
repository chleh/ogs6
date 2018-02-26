#include "CreateTCHSStokesMaterial.h"

namespace ProcessLib
{
namespace TCHSStokes
{
namespace Material
{
TCHSStokesMaterial createTCHSStokesMaterial(BaseLib::ConfigTree const& config)
{
    TCHSStokesMaterial material;

    material.fluid_density =
        createFluidDensity(config.getConfigSubtree("fluid_density"));
    material.fluid_viscosity =
        createFluidViscosity(config.getConfigSubtree("fluid_viscosity"));
    material.fluid_heat_capacity =
        createFluidHeatCapacity(config.getConfigSubtree("fluid_heat_capacity"));

    material.heat_conductivity =
        createHeatConductivity(config.getConfigSubtree("heat_conductivity"));

    material.diffusion_coefficient =
        ProcessLib::TES::createDiffusionCoefficient(
            config.getConfigSubtree("diffusion_coefficient"));

    material.fluid_momentum_production_coefficient =
        createFluidMomentumProductionCoefficient(
            config.getConfigSubtree("fluid_momentum_production_coefficient"));

    material.solid_heat_capacity =
        createSolidHeatCapacity(config.getConfigSubtree("solid_heat_capacity"));

    material.reaction_rate = MaterialLib::createReactionRate(
        config.getConfigSubtree("reaction_rate"));

    return material;
}

std::unordered_map<int, TCHSStokesMaterial> createTCHSStokesMaterials(
    BaseLib::ConfigTree const& config)
{
    std::unordered_map<int, TCHSStokesMaterial> materials;
    for (auto const c : config.getConfigSubtreeList("material"))
    {
        auto const id = c.getConfigAttribute<int>("id");
        auto it_succ = materials.emplace(std::piecewise_construct,
                                         std::forward_as_tuple(id),
                                         std::forward_as_tuple());
        if (!it_succ.second)
        {
            OGS_FATAL("Id %d already present.", id);
        }

        it_succ.first->second = createTCHSStokesMaterial(c);
    }
    return materials;
}

std::unique_ptr<FluidDensity> createFluidDensity(
    BaseLib::ConfigTree const& config)
{
    return nullptr;
}

std::unique_ptr<FluidViscosity> createFluidViscosity(
    BaseLib::ConfigTree const& config)
{
    return nullptr;
}

std::unique_ptr<FluidHeatCapacity> createFluidHeatCapacity(
    BaseLib::ConfigTree const& config)
{
    return nullptr;
}

std::unique_ptr<HeatConductivity> createHeatConductivity(
    BaseLib::ConfigTree const& config)
{
    return nullptr;
}

std::unique_ptr<FluidMomentumProductionCoefficient>
createFluidMomentumProductionCoefficient(BaseLib::ConfigTree const& config)
{
    return nullptr;
}

std::unique_ptr<SolidHeatCapacity> createSolidHeatCapacity(
    BaseLib::ConfigTree const& config)
{
    return nullptr;
}

}  // namespace Material

}  // namespace TCHSStokes

}  // namespace ProcessLib
