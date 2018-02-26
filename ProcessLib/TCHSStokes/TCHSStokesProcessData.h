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

#include "Material/TCHSStokesMaterial.h"
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
namespace TCHSStokes
{
template <int VelocityDim>
struct TCHSStokesProcessData
{
    enum
    {
        MATID_VOID = 0,
        MATID_BED = 1
    };

    TCHSStokesProcessData(
        MeshLib::PropertyVector<int> const& material_ids_,
        std::unordered_map<int, Material::TCHSStokesMaterial>&& materials_)
        : material_ids(material_ids_), materials(std::move(materials_))
    {
    }

    MeshLib::PropertyVector<int> const& material_ids;
    std::unordered_map<int, Material::TCHSStokesMaterial> materials;

    MeshLib::PropertyVector<double>* mesh_prop_nodal_p = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_nodal_T = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_nodal_xmV = nullptr;
};

}  // namespace TCHSStokes
}  // namespace ProcessLib
