/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MeshLib/PropertyVector.h"
#include "ProcessLib/Parameter/Parameter.h"

#include "ProcessLib/IncompressibleStokesBrinkman/Material/FluidViscosity.h"

#include <memory>
#include <utility>

#include <Eigen/Dense>

namespace MaterialLib
{
namespace Solids
{
template <int DisplacementDim>
struct MechanicsBase;
}
}
namespace ProcessLib
{
namespace TCHSStokes
{
template <int DisplacementDim>
struct TCHSStokesProcessData
{
    enum
    {
        MATID_VOID = 0,
        MATID_BED = 1
    };

    TCHSStokesProcessData(
        MeshLib::PropertyVector<int> const& material_ids_,
        double const pellet_diameter_,
        double const bed_radius_,
        double const average_darcy_velocity_,
        double const homogeneous_porosity_,
        ProcessLib::Parameter<double> const& fluid_density_,
        ProcessLib::Parameter<double> const& fluid_viscosity_,
        ProcessLib::Parameter<double> const& porosity_,
        std::unique_ptr<EffectiveFluidViscosity>&& effective_fluid_viscosity_)
        : material_ids(material_ids_),
          pellet_diameter(pellet_diameter_),
          bed_radius(bed_radius_),
          average_darcy_velocity(average_darcy_velocity_),
          homogeneous_porosity(homogeneous_porosity_),
          fluid_density(fluid_density_),
          fluid_viscosity(fluid_viscosity_),
          porosity(porosity_),
          effective_fluid_viscosity(std::move(effective_fluid_viscosity_))
    {
    }

#if 0
    TCHSStokesProcessData(
        TCHSStokesProcessData&& other)
        : material{std::move(other.material)},
          intrinsic_permeability(other.intrinsic_permeability),
          specific_storage(other.specific_storage),
          fluid_viscosity(other.fluid_viscosity),
          fluid_density(other.fluid_density),
          biot_coefficient(other.biot_coefficient),
          porosity(other.porosity),
          solid_density(other.solid_density),
          specific_body_force(other.specific_body_force),
          dt(other.dt),
          t(other.t)
    {
    }

    //! Copies are forbidden.
    TCHSStokesProcessData(
        TCHSStokesProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(TCHSStokesProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(TCHSStokesProcessData&&) = delete;
#endif

    MeshLib::PropertyVector<int> const& material_ids;

    double const pellet_diameter;
    double const bed_radius;
    double const average_darcy_velocity;
    double const homogeneous_porosity;
    ProcessLib::Parameter<double> const& fluid_density;
    ProcessLib::Parameter<double> const& fluid_viscosity;
    ProcessLib::Parameter<double> const& porosity;

    std::unique_ptr<EffectiveFluidViscosity> effective_fluid_viscosity;

    double dt = 0.0;
    double t = 0.0;

    MeshLib::PropertyVector<double>* mesh_prop_nodal_p = nullptr;
};

}  // namespace TCHSStokes
}  // namespace ProcessLib
