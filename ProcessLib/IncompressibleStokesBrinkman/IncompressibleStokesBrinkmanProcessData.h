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
    IncompressibleStokesBrinkmanProcessData(
        Parameter<double> const& porosity_,
        Parameter<double> const& mu_eff_,
        Parameter<double> const& lambda_eff_,
        Parameter<double> const& f_1_,
        Parameter<double> const& f_2_)
        : porosity(porosity_),
          mu_eff(mu_eff_),
          lambda_eff(lambda_eff_),
          f_1(f_1_),
          f_2(f_2_)
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

    Parameter<double> const& porosity;
    Parameter<double> const& mu_eff;
    Parameter<double> const& lambda_eff;
    Parameter<double> const& f_1;
    Parameter<double> const& f_2;
    double dt = 0.0;
    double t = 0.0;

    MeshLib::PropertyVector<double>* mesh_prop_nodal_p = nullptr;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace IncompressibleStokesBrinkman
}  // namespace ProcessLib
