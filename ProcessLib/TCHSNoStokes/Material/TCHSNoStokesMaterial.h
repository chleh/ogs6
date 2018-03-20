#pragma once

#include <cassert>
#include <memory>

#include "MaterialLib/ReactionKinetics/ReactionRate.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "ProcessLib/IncompressibleStokesBrinkman/Material/FluidViscosity.h"
#include "ProcessLib/TCHSStokes/Material/TCHSStokesMaterial.h"
#include "ProcessLib/TES/Material/DiffusionCoefficient.h"
#include "ProcessLib/TES/TESOGS5MaterialModels.h"

namespace ProcessLib
{
namespace TCHSNoStokes
{
namespace Material
{
struct Permeability
{
    virtual double getPermeability(double r) const = 0;
    virtual ~Permeability() = default;
};

struct ConstantPermeability final : Permeability
{
    explicit ConstantPermeability(double value) : _value(value) {}

    double getPermeability(double /*r*/) const override { return _value; }

private:
    double const _value;
};

struct RadialPermeability final : Permeability
{
    explicit RadialPermeability(
        std::unique_ptr<MathLib::PiecewiseLinearInterpolation>&& radial_profile)
        : _radial_profile(std::move(radial_profile))
    {
    }

    double getPermeability(double r) const override
    {
        return _radial_profile->getValue(r);
    }

private:
    std::unique_ptr<MathLib::PiecewiseLinearInterpolation> const
        _radial_profile;
};

struct TCHSNoStokesMaterial
{
    std::unique_ptr<ProcessLib::TCHSStokes::Material::FluidDensity>
        fluid_density;
    std::unique_ptr<ProcessLib::TCHSStokes::Material::FluidViscosity>
        fluid_viscosity;
    std::unique_ptr<ProcessLib::TCHSStokes::Material::FluidHeatCapacity>
        fluid_heat_capacity;
    std::unique_ptr<ProcessLib::TCHSStokes::Material::HeatConductivity>
        heat_conductivity;
    std::unique_ptr<ProcessLib::TCHSStokes::Material::MassDispersion>
        mass_dispersion;
    std::unique_ptr<ProcessLib::TCHSStokes::Material::SolidHeatCapacity>
        solid_heat_capacity;
    std::unique_ptr<MaterialLib::ReactiveSolidModel> reactive_solid;
    std::unique_ptr<MaterialLib::ReactionRate> reaction_rate;
    std::unique_ptr<ProcessLib::TCHSStokes::Material::Porosity> porosity;

    std::unique_ptr<ProcessLib::ReynoldsNumber> reynolds_number;
    std::unique_ptr<ProcessLib::TCHSStokes::Material::PecletNumberHeat>
        peclet_number_heat;

    std::unique_ptr<Permeability> permeability;

    double molar_mass_reactive = std::numeric_limits<double>::quiet_NaN();
    double molar_mass_inert = std::numeric_limits<double>::quiet_NaN();
};

}  // namespace Material

}  // namespace TCHSNoStokes

}  // namespace ProcessLib
