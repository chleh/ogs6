/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <unordered_map>
#include <utility>

#include <Eigen/Dense>

#include "MeshLib/PropertyVector.h"
#include "ProcessLib/Parameter/Parameter.h"

#include "Material/TCHSNoStokesMaterial.h"
#include "MaterialLib/ReactionKinetics/ReactionRate.h"
#include "ProcessLib/IncompressibleStokesBrinkman/Material/FluidViscosity.h"

namespace MaterialLib
{
namespace Solids
{
template <int VelocityDim>
struct MechanicsBase;
}
}
namespace ProcessLib
{
namespace TCHSNoStokes
{
template <int VelocityDim>
struct TCHSNoStokesProcessData
{
    TCHSNoStokesProcessData(
        MeshLib::PropertyVector<int> const& material_ids_,
        std::unordered_map<int, Material::TCHSNoStokesMaterial>&& materials_,
        std::size_t velocity_probe_node_id_,
        double darcy_velocity_center_)
        : material_ids(material_ids_),
          materials(std::move(materials_)),
          velocity_probe_node_id(velocity_probe_node_id_),
          darcy_velocity_center(darcy_velocity_center_)
    {
    }

    MeshLib::PropertyVector<int> const& material_ids;
    std::unordered_map<int, Material::TCHSNoStokesMaterial> materials;

    double delta_t = std::numeric_limits<double>::quiet_NaN();

    // for heat conductivity computation
    std::size_t const velocity_probe_node_id;
    double probed_pressure = std::numeric_limits<double>::quiet_NaN();
    double probed_temperature = std::numeric_limits<double>::quiet_NaN();
    const double darcy_velocity_center;

    MeshLib::PropertyVector<double>* mesh_prop_cell_hat_rho_SR = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_cell_reaction_enthalpy = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_cell_rho_SR = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_cell_rho_GR = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_cell_lambda = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_cell_cpS = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_cell_cpG = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_cell_v_Darcy = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_cell_mass_flux = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_cell_vapour_mass_flux = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_cell_solid_mass = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_cell_heating_rate = nullptr;
};

}  // namespace TCHSNoStokes
}  // namespace ProcessLib
