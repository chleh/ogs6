/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_BOUNDARY_CONDITION_H_
#define PROCESS_LIB_BOUNDARY_CONDITION_H_

#include "BaseLib/ConfigTree.h"
#include "MeshGeoToolsLib/MeshNodeSearcher.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/NumericsConfig.h" // for GlobalIndexType

#include "DirichletBc.h"

namespace GeoLib
{
class GeoObject;
}

namespace ProcessLib
{
/// The UniformDirichletBoundaryCondition class describes a constant in space
/// and time Dirichlet boundary condition.
/// The expected parameter in the passed configuration is "value" which, when
/// not present defaults to zero.
class UniformDirichletBoundaryCondition
{
public:
    UniformDirichletBoundaryCondition(GeoLib::GeoObject const& geometry,
                                      BaseLib::ConfigTree const& config);

    UniformDirichletBoundaryCondition(GeoLib::GeoObject const& geometry,
                                      double value);

    std::unique_ptr<DirichletBoundaryCondition> getDirichletBoundaryCondition(
        MeshGeoToolsLib::MeshNodeSearcher& searcher,
        NumLib::LocalToGlobalIndexMap const& dof_table,
        int const variable_id,
        int const component_id);

private:
    double _value;
    GeoLib::GeoObject const& _geometry;
};

}  // namespace ProcessLib

#endif  // PROCESS_LIB_BOUNDARY_CONDITION_H_
