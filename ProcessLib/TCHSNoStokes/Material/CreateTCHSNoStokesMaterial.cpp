#include "CreateTCHSNoStokesMaterial.h"
#include "ProcessLib/TCHSStokes/Material/CreateTCHSStokesMaterial.h"

namespace ProcessLib
{
namespace TCHSNoStokes
{
namespace Material
{
TCHSNoStokesMaterial createTCHSNoStokesMaterial(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters)
{
    using namespace ProcessLib::TCHSStokes::Material;

    TCHSNoStokesMaterial material;

    material.fluid_density =
        createFluidDensity(config.getConfigSubtree("fluid_density"));
    material.fluid_viscosity =
        createFluidViscosity(config.getConfigSubtree("fluid_viscosity"));
    material.fluid_heat_capacity =
        createFluidHeatCapacity(config.getConfigSubtree("fluid_heat_capacity"));

    material.heat_conductivity =
        createHeatConductivity(config.getConfigSubtree("heat_conductivity"));

    material.mass_dispersion =
        createMassDispersion(config.getConfigSubtree("mass_dispersion"));

    material.solid_heat_capacity =
        createSolidHeatCapacity(config.getConfigSubtree("solid_heat_capacity"));

    material.reaction_rate = MaterialLib::createReactionRate(
        config.getConfigSubtree("reaction_rate"));

    material.reactive_solid = MaterialLib::createReactiveSolidModel(
        config.getConfigSubtree("reactive_solid"), parameters);

    material.porosity = createPorosity(config.getConfigSubtree("porosity"));

    material.molar_mass_reactive =
        config.getConfigParameter<double>("molar_mass_reactive");
    material.molar_mass_inert =
        config.getConfigParameter<double>("molar_mass_inert");

    material.reynolds_number = ProcessLib::createReynoldsNumber(
        config.getConfigSubtree("reynolds_number"));

    return material;
}

std::unordered_map<int, TCHSNoStokesMaterial> createTCHSNoStokesMaterials(
    BaseLib::ConfigTree const& config,
    const std::vector<std::unique_ptr<ParameterBase>>& parameters)
{
    std::unordered_map<int, TCHSNoStokesMaterial> materials;
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

        it_succ.first->second = createTCHSNoStokesMaterial(c, parameters);
    }
    return materials;
}

}  // namespace Material

}  // namespace TCHSNoStokes

}  // namespace ProcessLib
