#pragma once

#include <memory>

#include "MaterialLib/ReactionKinetics/ReactionRate.h"
#include "ProcessLib/TES/Material/DiffusionCoefficient.h"

namespace ProcessLib
{
namespace TCHSStokes
{
namespace Material
{
class FluidDensity
{
public:
    virtual ~FluidDensity() = default;
};

class FluidDensityMixtureWaterNitrogen final : public FluidDensity
{
};

class FluidViscosity
{
public:
    virtual ~FluidViscosity() = default;
};

class FluidViscosityMixtureWaterNitrogen final : public FluidViscosity
{
};

class FluidHeatCapacity
{
public:
    virtual ~FluidHeatCapacity() = default;
};

class FluidHeatCapacityMixtureWaterNitrogen final : public FluidHeatCapacity
{
};

class HeatConductivity
{
public:
    virtual ~HeatConductivity() = default;
};

class HeatConductivityLambdaR final : public HeatConductivity
{
};

class FluidMomentumProductionCoefficient
{
public:
    virtual ~FluidMomentumProductionCoefficient() = default;
};

class FluidMomentumProductionCoefficientErgun final
    : public FluidMomentumProductionCoefficient
{
};

class SolidHeatCapacity
{
public:
    virtual ~SolidHeatCapacity() = default;
};

class SolidHeatCapacityConstant final : public SolidHeatCapacity
{
};

struct TCHSStokesMaterial
{
    std::unique_ptr<FluidDensity> fluid_density;
    std::unique_ptr<FluidViscosity> fluid_viscosity;
    std::unique_ptr<FluidHeatCapacity> fluid_heat_capacity;
    std::unique_ptr<HeatConductivity> heat_conductivity;
    std::unique_ptr<ProcessLib::TES::DiffusionCoefficient>
        diffusion_coefficient;
    std::unique_ptr<FluidMomentumProductionCoefficient>
        fluid_momentum_production_coefficient;
    std::unique_ptr<SolidHeatCapacity> solid_heat_capacity;
    std::unique_ptr<MaterialLib::ReactionRate> reaction_rate;
};

}  // namespace Material

}  // namespace TCHSStokes

}  // namespace ProcessLib
