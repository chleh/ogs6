/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <utility>

#include <Eigen/Eigen>

#if HAVE_UG
#pragma message "HAVE_UG set ##### SmallDefPcsData"
#else
#pragma message "HAVE_UG not set ##### SmallDefPcsData"
#endif

#include "MeshLib/DUNEMesh.h"
#include "ProcessLib/Parameter/Parameter.h"

namespace MaterialLib
{
namespace Solids
{
template <int DisplacementDim>
struct MechanicsBase;
}
}  // namespace MaterialLib
namespace ProcessLib
{
namespace SmallDeformationDUNE
{
template <int DisplacementDim>
struct SmallDeformationProcessData
{
    SmallDeformationProcessData(
        std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>&&
            material,
        Parameter<double> const& solid_density_,
        Eigen::Matrix<double, DisplacementDim, 1>
            specific_body_force_,
        double const reference_temperature_,
        MeshLib::DUNEMesh<DisplacementDim>& grid_)
        : material{std::move(material)},
          solid_density(solid_density_),
          specific_body_force(std::move(specific_body_force_)),
          reference_temperature(reference_temperature_),
          grid(grid_)
    {
    }

    SmallDeformationProcessData(SmallDeformationProcessData&& other)
        : material{std::move(other.material)},
          solid_density(other.solid_density),
          specific_body_force(other.specific_body_force),
          dt{other.dt},
          t{other.t},
          reference_temperature(other.reference_temperature),
          grid{other.grid}
    {
    }

    //! Copies are forbidden.
    SmallDeformationProcessData(SmallDeformationProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(SmallDeformationProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(SmallDeformationProcessData&&) = delete;

    std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>
        material;
    /// Solid's density. A scalar quantity, Parameter<double>.
    Parameter<double> const& solid_density;
    /// Specific body forces applied to the solid.
    /// It is usually used to apply gravitational forces.
    /// A vector of displacement dimension's length.
    Eigen::Matrix<double, DisplacementDim, 1> const specific_body_force;
    double dt = 0;
    double t = 0;
    double const reference_temperature;

    MeshLib::DUNEMesh<DisplacementDim>& grid;  // maybe gridView?
};

}  // namespace SmallDeformationDUNE
}  // namespace ProcessLib
