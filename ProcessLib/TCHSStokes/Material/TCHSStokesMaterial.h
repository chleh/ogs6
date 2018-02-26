#pragma once

#include <memory>

#include "MaterialLib/ReactionKinetics/ReactionRate.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "ProcessLib/IncompressibleStokesBrinkman/Material/FluidViscosity.h"
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
    double operator()()
    {
        // TODO remove
        return 0.0;
    }
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
    virtual double getPorosity(double r) = 0;
    virtual double getDPorosityDr(double r) = 0;
    ~Porosity() = default;
};

class PorosityConstant : public Porosity
{
public:
    PorosityConstant(double value) : _value(value) {}

    double getPorosity(double /*r*/) override { return _value; }
    double getDPorosityDr(double /*r*/) override { return 0.0; }

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

    double getPorosity(double r) override
    {
        auto const exp_term =
            std::exp(-5.0 * (_bed_radius - r) / _pellet_diameter);
        return _homogeneous_porosity + _homogeneous_porosity * 1.36 * exp_term;
    }

    double getDPorosityDr(double r) override
    {
        auto const exp_term =
            std::exp(-5.0 * (_bed_radius - r) / _pellet_diameter);
        return _homogeneous_porosity * 1.36 * 5.0 / _pellet_diameter * exp_term;
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
    std::unique_ptr<EffectiveFluidViscosity> effective_fluid_viscosity;
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
