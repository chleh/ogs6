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
struct TCHSNoStokesMaterial
{
    std::unique_ptr<ProcessLib::TCHSStokes::Material::FluidDensity>
        fluid_density;
    std::unique_ptr<ProcessLib::TCHSStokes::Material::FluidViscosity>
        fluid_viscosity;
    std::unique_ptr<EffectiveFluidViscosity> effective_fluid_viscosity;
    std::unique_ptr<ProcessLib::TCHSStokes::Material::FluidHeatCapacity>
        fluid_heat_capacity;
    std::unique_ptr<ProcessLib::TCHSStokes::Material::HeatConductivity>
        heat_conductivity;
    std::unique_ptr<ProcessLib::TCHSStokes::Material::MassDispersion>
        mass_dispersion;
    std::unique_ptr<
        ProcessLib::TCHSStokes::Material::FluidMomentumProductionCoefficient>
        fluid_momentum_production_coefficient;
    std::unique_ptr<ProcessLib::TCHSStokes::Material::SolidHeatCapacity>
        solid_heat_capacity;
    std::unique_ptr<MaterialLib::ReactiveSolidModel> reactive_solid;
    std::unique_ptr<MaterialLib::ReactionRate> reaction_rate;
    std::unique_ptr<ProcessLib::TCHSStokes::Material::Porosity> porosity;

    std::unique_ptr<ProcessLib::ReynoldsNumber> reynolds_number;
    std::unique_ptr<ProcessLib::TCHSStokes::Material::PecletNumberHeat>
        peclet_number_heat;

    double molar_mass_reactive = std::numeric_limits<double>::quiet_NaN();
    double molar_mass_inert = std::numeric_limits<double>::quiet_NaN();
};

}  // namespace Material

}  // namespace TCHSNoStokes

}  // namespace ProcessLib
