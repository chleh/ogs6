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
    virtual double getHeatCapacity(double rho_SR, double T) const = 0;
    virtual ~SolidHeatCapacity() = default;
};

class SolidHeatCapacityZeoliteAWaterVucelic final : public SolidHeatCapacity
{
public:
    SolidHeatCapacityZeoliteAWaterVucelic(double cp_zeo_dry, double rho_SR_dry)
        : _cp_zeo_dry(cp_zeo_dry), _rho_SR_dry(rho_SR_dry)
    {
    }

    double getHeatCapacity(double rho_SR, double T) const override
    {
        double const cp_min_water = 836.0;  // J/kg/K
        double const A = 3.762e-2;
        double const B = 3.976e-4;
        double const sqrtB = std::sqrt(B);
        double const T_0 = 335;    // K
        double const T_min = 220;  // K; TODO check, estimated from Vučelić, V.,
                                   // Vučelić, D., 1983. The heat capacity of
                                   // water near solid surfaces. Chemical
                                   // Physics Letters 102, 371–374.
                                   // doi:10.1016/0009-2614(83)87058-4
        double const sqrtpi = std::sqrt(boost::math::constants::pi<double>());

        double const loading = rho_SR / _rho_SR_dry - 1.0;

        double const cp_water =
            cp_min_water +
            0.5 * sqrtpi * A / sqrtB *
                (std::erf(sqrtB * (T - T_0)) - std::erf(sqrtB * (T_min - T_0)));

        return 1.0 / (1.0 + loading) * _cp_zeo_dry +
               loading / (1.0 + loading) * cp_water;
    }

private:
    double const _cp_zeo_dry;
    double const _rho_SR_dry;
};

class SolidHeatCapacityConstant final : public SolidHeatCapacity
{
public:
    SolidHeatCapacityConstant(double value) : _value(value) {}

    double getHeatCapacity(double /*rho_SR*/, double /*T*/) const override
    {
        return _value;
    }

private:
    double const _value;
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

    double molar_mass_reactive = std::numeric_limits<double>::quiet_NaN();
    double molar_mass_inert = std::numeric_limits<double>::quiet_NaN();
};

}  // namespace Material

}  // namespace TCHSStokes

}  // namespace ProcessLib
