#pragma once

#include <memory>

#include "MaterialLib/ReactionKinetics/ReactionRate.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "ProcessLib/IncompressibleStokesBrinkman/Material/FluidViscosity.h"
#include "ProcessLib/TES/Material/DiffusionCoefficient.h"
#include "ProcessLib/TES/TESOGS5MaterialModels.h"

namespace ProcessLib
{
namespace TCHSStokes
{
namespace Material
{
class FluidDensity
{
public:
    virtual double getDensity(const double p,
                              const double T,
                              const double x_mV) const = 0;

    virtual ~FluidDensity() = default;
};

class FluidDensityMixtureWaterNitrogen final : public FluidDensity
{
public:
    double getDensity(const double p,
                      const double T,
                      const double x_mV) const override
    {
        return ProcessLib::TES::fluid_density(p, T, x_mV);
    }
};

class FluidViscosity
{
public:
    virtual double getViscosity(const double p,
                                const double T,
                                const double x) const = 0;
    virtual ~FluidViscosity() = default;
};

class FluidViscosityMixtureWaterNitrogen final : public FluidViscosity
{
public:
    double getViscosity(const double p,
                        const double T,
                        const double x) const override
    {
        return ProcessLib::TES::fluid_viscosity(p, T, x);
    }
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
    virtual double getCoeffOfV(double porosity, double viscosity) const = 0;
    virtual double getCoeffOfVSquared(double porosity,
                                      double fluid_density) const = 0;

    virtual ~FluidMomentumProductionCoefficient() = default;
};

class FluidMomentumProductionCoefficientErgun final
    : public FluidMomentumProductionCoefficient
{
public:
    FluidMomentumProductionCoefficientErgun(const double pellet_diameter)
        : _pellet_diameter(pellet_diameter)
    {
    }

    double getCoeffOfV(double porosity, double viscosity) const override
    {
        auto const poro3 = boost::math::pow<3>(porosity);
        return 150.0 * boost::math::pow<2>(1.0 - porosity) / poro3 * viscosity /
               _pellet_diameter / _pellet_diameter;
    }
    double getCoeffOfVSquared(double porosity,
                              double fluid_density) const override
    {
        auto const poro3 = boost::math::pow<3>(porosity);
        return 1.75 * (1.0 - porosity) / poro3 * fluid_density /
               _pellet_diameter;
    }

private:
    double const _pellet_diameter;
};

class FluidMomentumProductionCoefficientZero final
    : public FluidMomentumProductionCoefficient
{
    double getCoeffOfV(double /*porosity*/, double /*viscosity*/) const override
    {
        return 0;
    }
    double getCoeffOfVSquared(double /*porosity*/,
                              double /*fluid_density*/) const override
    {
        return 0;
    }
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
    virtual double getPorosity(double r) const = 0;
    virtual double getDPorosityDr(double r) const = 0;
    ~Porosity() = default;
};

class PorosityConstant : public Porosity
{
public:
    PorosityConstant(double value) : _value(value) {}

    double getPorosity(double /*r*/) const override { return _value; }
    double getDPorosityDr(double /*r*/) const override { return 0.0; }

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

    double getPorosity(double r) const override
    {
        auto const exp_term =
            std::exp(-5.0 * (_bed_radius - r) / _pellet_diameter);
        return _homogeneous_porosity + _homogeneous_porosity * 1.36 * exp_term;
    }

    double getDPorosityDr(double r) const override
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
    std::unique_ptr<MaterialLib::ReactiveSolidModel> reactive_solid;
    std::unique_ptr<MaterialLib::ReactionRate> reaction_rate;
    std::unique_ptr<Porosity> porosity;
};

}  // namespace Material

}  // namespace TCHSStokes

}  // namespace ProcessLib
