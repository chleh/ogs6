#include "CreateTCHSStokesMaterial.h"

namespace ProcessLib
{
namespace TCHSStokes
{
namespace Material
{
TCHSStokesMaterial createTCHSStokesMaterial(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters)
{
    TCHSStokesMaterial material;

    material.fluid_density =
        createFluidDensity(config.getConfigSubtree("fluid_density"));
    material.fluid_viscosity =
        createFluidViscosity(config.getConfigSubtree("fluid_viscosity"));
    material.effective_fluid_viscosity =
        ProcessLib::createEffectiveFluidViscosity(
            config.getConfigSubtree("effective_fluid_viscosity"));
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

    material.reactive_solid = MaterialLib::createReactiveSolidModel(
        config.getConfigSubtree("reactive_solid"), parameters);

    material.porosity = createPorosity(config.getConfigSubtree("porosity"));

    material.molar_mass_reactive =
        config.getConfigParameter<double>("molar_mass_reactive");
    material.molar_mass_inert =
        config.getConfigParameter<double>("molar_mass_inert");

    return material;
}

std::unordered_map<int, TCHSStokesMaterial> createTCHSStokesMaterials(
    BaseLib::ConfigTree const& config,
    const std::vector<std::unique_ptr<ParameterBase>>& parameters)
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

        it_succ.first->second = createTCHSStokesMaterial(c, parameters);
    }
    return materials;
}

std::unique_ptr<FluidDensity> createFluidDensity(
    BaseLib::ConfigTree const& config)
{
    auto const type = config.getConfigParameter<std::string>("type");
    if (type == "MixtureWaterNitrogen")
        return std::make_unique<FluidDensityMixtureWaterNitrogen>();

    OGS_FATAL("Unknown fluid density model `%s'.", type.c_str());
}

std::unique_ptr<FluidViscosity> createFluidViscosity(
    BaseLib::ConfigTree const& config)
{
    auto const type = config.getConfigParameter<std::string>("type");
    if (type == "MixtureWaterNitrogen")
        return std::make_unique<FluidViscosityMixtureWaterNitrogen>();

    OGS_FATAL("Unknown fluid viscosity model `%s'.", type.c_str());
}

std::unique_ptr<FluidHeatCapacity> createFluidHeatCapacity(
    BaseLib::ConfigTree const& config)
{
    auto const type = config.getConfigParameter<std::string>("type");
    if (type == "MixtureWaterNitrogen")
        return std::make_unique<FluidHeatCapacityMixtureWaterNitrogen>();

    OGS_FATAL("Unknown fluid heat capacity model `%s'.", type.c_str());
}

std::unique_ptr<HeatConductivity> createHeatConductivity(
    BaseLib::ConfigTree const& config)
{
    auto const type = config.getConfigParameter<std::string>("type");
    if (type == "LambdaRModelNoRadiation")
    {
        auto lambda_f = createHeatConductivity(
            config.getConfigSubtree("fluid_heat_conductivity"));
        auto lambda_p = createHeatConductivity(
            config.getConfigSubtree("solid_heat_conductivity"));

        auto const pellet_diameter =
            config.getConfigParameter<double>("pellet_diameter");
        auto const bed_radius = config.getConfigParameter<double>("bed_radius");

        return std::make_unique<HeatConductivityLambdaRNoRadiation>(
            bed_radius,
            pellet_diameter,
            std::move(lambda_f),
            std::move(lambda_p));
    }
    if (type == "MixtureWaterNitrogen")
        return std::make_unique<HeatConductivityMixtureWaterNitrogen>();
    if (type == "Constant")
    {
        auto const value = config.getConfigParameter<double>("value");
        return std::make_unique<HeatConductivityConstant>(value);
    }

    OGS_FATAL("Unknown heat conductivity model `%s'.", type.c_str());
}

std::unique_ptr<FluidMomentumProductionCoefficient>
createFluidMomentumProductionCoefficient(BaseLib::ConfigTree const& config)
{
    auto const type = config.getConfigParameter<std::string>("type");
    if (type == "Ergun")
    {
        auto const pellet_diameter =
            config.getConfigParameter<double>("pellet_diameter");

        return std::make_unique<FluidMomentumProductionCoefficientErgun>(
            pellet_diameter);
    }
    if (type == "Zero")
        return std::make_unique<FluidMomentumProductionCoefficientZero>();

    OGS_FATAL("Unknown fluid momentum production model `%s'.", type.c_str());
}

std::unique_ptr<SolidHeatCapacity> createSolidHeatCapacity(
    BaseLib::ConfigTree const& config)
{
    auto const type = config.getConfigParameter<std::string>("type");
    if (type == "Constant")
    {
        auto const value = config.getConfigParameter<double>("value");
        return std::make_unique<SolidHeatCapacityConstant>(value);
    }
    if (type == "ZeoliteAWaterVucelic")
    {
        auto const rho_SR_dry =
            config.getConfigParameter<double>("adsorbent_density_dry");
        auto const cp_zeo_dry =
            config.getConfigParameter<double>("adsorbent_heat_capacity_dry");
        return std::make_unique<SolidHeatCapacityZeoliteAWaterVucelic>(
            rho_SR_dry, cp_zeo_dry);
    }

    OGS_FATAL("Unknown solid heat capacity model `%s'.", type.c_str());
}

std::unique_ptr<Porosity> createPorosity(BaseLib::ConfigTree const& config)
{
    auto const type = config.getConfigParameter<std::string>("type");
    if (type == "Constant")
    {
        auto const value = config.getConfigParameter<double>("value");
        return std::make_unique<PorosityConstant>(value);
    }
    if (type == "RadialProfileGiese")
    {
        auto const pellet_diameter =
            config.getConfigParameter<double>("pellet_diameter");
        auto const bed_radius = config.getConfigParameter<double>("bed_radius");
        auto const homogeneous_porosity =
            config.getConfigParameter<double>("homogeneous_porosity");
        return std::make_unique<PorosityRadialProfileGiese>(
            bed_radius, pellet_diameter, homogeneous_porosity);
    }

    OGS_FATAL("Unknown porosity model `%s'.", type.c_str());
}

}  // namespace Material

}  // namespace TCHSStokes

}  // namespace ProcessLib
