#pragma once

#include <memory>

#include "MaterialLib/ReactionKinetics/ReactionRate.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
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
    double operator()()
    {
        // TODO remove
        return 0.0;
    }
    virtual ~FluidViscosity() = default;
};

class FluidViscosityMixtureWaterNitrogen final : public FluidViscosity
{
};

class FluidViscosityEffectiveGiese final : public FluidViscosity
{
public:
    FluidViscosityEffectiveGiese(
        std::unique_ptr<FluidViscosity>&& model,
        MathLib::PiecewiseLinearInterpolation&& average_darcy_velocity,
        double pellet_diameter)
        : _model(std::move(model)),
          _average_darcy_velocity(std::move(average_darcy_velocity)),
          _pellet_diameter(pellet_diameter)
    {
    }

private:
    std::unique_ptr<FluidViscosity> _model;
    MathLib::PiecewiseLinearInterpolation _average_darcy_velocity;
    const double _pellet_diameter;
};

class FluidHeatCapacity
{
public:
    double operator()()
    {
        // TODO remove
        return 0.0;
    }
    virtual ~FluidHeatCapacity() = default;
};

class FluidHeatCapacityMixtureWaterNitrogen final : public FluidHeatCapacity
{
};

class HeatConductivity
{
public:
    double operator()()
    {
        // TODO remove
        return 0.0;
    }
    virtual ~HeatConductivity() = default;
};

class HeatConductivityLambdaR final : public HeatConductivity
{
};

class HeatConductivityMixtureWaterNitrogen final : public HeatConductivity
{
};

class FluidMomentumProductionCoefficient
{
public:
    double operator()()
    {
        // TODO remove
        return 0.0;
    }
    virtual ~FluidMomentumProductionCoefficient() = default;
};

class FluidMomentumProductionCoefficientErgun final
    : public FluidMomentumProductionCoefficient
{
};

class FluidMomentumProductionCoefficientZero final
    : public FluidMomentumProductionCoefficient
{
};

class SolidHeatCapacity
{
public:
    double operator()()
    {
        // TODO remove
        return 0.0;
    }
    virtual ~SolidHeatCapacity() = default;
};

class SolidHeatCapacityConstant final : public SolidHeatCapacity
{
public:
    SolidHeatCapacityConstant(double value) : _value(value) {}

private:
    double const _value;
};

class SolidHeatCapacityDensityDependent final : public SolidHeatCapacity
{
};

class Porosity
{
public:
    double operator()()
    {
        // TODO remove
        return 0.0;
    }
    ~Porosity() = default;
};

class PorosityConstant : public Porosity
{
public:
    PorosityConstant(double value) : _value(value) {}

private:
    double const _value;
};

class PorosityRadialProfileGiese : public Porosity
{
public:
    PorosityRadialProfileGiese(const double bed_radius,
                               const double pellet_diameter,
                               const double homogeneous_porosity)
        : _bed_radius(bed_radius),
          _pellet_diameter(pellet_diameter),
          _homogeneous_porosity(homogeneous_porosity)
    {
    }

private:
    const double _bed_radius;
    const double _pellet_diameter;
    const double _homogeneous_porosity;
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
    std::unique_ptr<Porosity> porosity;
};

}  // namespace Material

}  // namespace TCHSStokes

}  // namespace ProcessLib
