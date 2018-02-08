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
namespace IncompressibleStokesBrinkman
{
template <int DisplacementDim>
struct IncompressibleStokesBrinkmanProcessData
{
    enum
    {
        MATID_VOID = 0,
        MATID_BED = 1
    };

    IncompressibleStokesBrinkmanProcessData(
        Parameter<int> const& materialIDs_,
        double const pellet_diameter_,
        double const bed_radius_,
        double const average_darcy_velocity_,
        double const fluid_density_,
        double const fluid_viscosity_)
        : materialIDs(materialIDs_),
          pellet_diameter(pellet_diameter_),
          bed_radius(bed_radius_),
          average_darcy_velocity(average_darcy_velocity_),
          fluid_density(fluid_density_),
          fluid_viscosity(fluid_viscosity_)
    {
    }

#if 0
    IncompressibleStokesBrinkmanProcessData(
        IncompressibleStokesBrinkmanProcessData&& other)
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
    IncompressibleStokesBrinkmanProcessData(
        IncompressibleStokesBrinkmanProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(IncompressibleStokesBrinkmanProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(IncompressibleStokesBrinkmanProcessData&&) = delete;
#endif

    Parameter<int> const& materialIDs;

    double const pellet_diameter;
    double const bed_radius;
    double const average_darcy_velocity;
    double const fluid_density;
    double const fluid_viscosity;

    double dt = 0.0;
    double t = 0.0;

    MeshLib::PropertyVector<double>* mesh_prop_nodal_p = nullptr;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace IncompressibleStokesBrinkman
}  // namespace ProcessLib
