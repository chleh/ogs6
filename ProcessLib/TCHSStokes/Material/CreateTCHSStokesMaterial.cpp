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

    material.porosity = createPorosity(config.getConfigSubtree("porosity"));

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
    if (type == "EffectiveGiese")
    {
        auto model =
            createFluidViscosity(config.getConfigSubtree("fluid_viscosity"));

        auto const c = config.getConfigSubtree("average_darcy_velocity");
        MathLib::PiecewiseLinearInterpolation average_darcy_velocity{
            c.getConfigParameter<std::vector<double>>("times"),
            c.getConfigParameter<std::vector<double>>("values")};

        auto const pellet_diameter =
            config.getConfigParameter<double>("pellet_diameter");

        return std::make_unique<FluidViscosityEffectiveGiese>(
            std::move(model),
            std::move(average_darcy_velocity),
            pellet_diameter);
    }

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
    if (type == "LambdaRModel")
    {
        config.ignoreConfigParameter("fluid_heat_conductivity");
        config.ignoreConfigParameter("solid_heat_conductivity");
        return std::make_unique<HeatConductivityLambdaR>();
    }
    if (type == "MixtureWaterNitrogen")
        return std::make_unique<HeatConductivityMixtureWaterNitrogen>();

    OGS_FATAL("Unknown heat conductivity model `%s'.", type.c_str());
}

std::unique_ptr<FluidMomentumProductionCoefficient>
createFluidMomentumProductionCoefficient(BaseLib::ConfigTree const& config)
{
    auto const type = config.getConfigParameter<std::string>("type");
    if (type == "Ergun")
        return std::make_unique<FluidMomentumProductionCoefficientErgun>();
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
    if (type == "DensityDependent")
    {
        return std::make_unique<SolidHeatCapacityDensityDependent>();
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
