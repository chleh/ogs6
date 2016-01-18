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

#include <algorithm>
#include <vector>

#include "logog/include/logog.hpp"

#include "NumericsConfig.h" // for GlobalIndexType

#include "BaseLib/ConfigTreeNew.h"
#include "AssemblerLib/LocalToGlobalIndexMap.h"
#include "MeshGeoToolsLib/MeshNodeSearcher.h"

#include "ProcessLib/VariableTransformation.h"

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
    UniformDirichletBoundaryCondition(GeoLib::GeoObject const* const geometry,
                                      BaseLib::ConfigTreeNew const& config)
        : _geometry(geometry)
    {
        DBUG("Constructing UniformDirichletBoundaryCondition from config.");
        config.checkConfParam("type", "UniformDirichlet");

        _value = config.getConfParam<double>("value");
        DBUG("Using value %g", _value);
    }

    /// Initialize Dirichlet type boundary conditions.
    /// Fills in global_ids of the particular geometry of the boundary condition
    /// and the corresponding values.
    /// The ids are appended to the global_ids and the values are filled with
    /// the constant _value.
    void initialize(
            MeshGeoToolsLib::MeshNodeSearcher& searcher,
            AssemblerLib::LocalToGlobalIndexMap const& dof_table,
            std::size_t component_id,
            std::vector<GlobalIndexType>& global_ids,
            std::vector<double>& values)
    {
        // Find nodes' ids on the given mesh on which this boundary condition
        // is defined.
        std::vector<std::size_t> ids = searcher.getMeshNodeIDs(*_geometry);

        // convert mesh node ids to global index for the given component
        global_ids.reserve(global_ids.size() + ids.size());
        values.reserve(values.size() + ids.size());
        for (auto& id : ids)
        {
            MeshLib::Location l(searcher.getMeshId(),
                                MeshLib::MeshItemType::Node,
                                id);
            // TODO: that might be slow, but only done once
            const auto g_idx = dof_table.getGlobalIndex(l, component_id);
            // For the DDC approach (e.g. with PETSc option), the negative
            // index of g_idx means that the entry by that index is a ghost one,
            // which should be dropped. Especially for PETSc routines MatZeroRows
            // and MatZeroRowsColumns, which are called to apply the Dirichlet BC,
            // the negative index is not accepted like other matrix or vector
            // PETSc routines. Therefore, the following if-condition is applied.
            if (g_idx >= 0)
            {
                global_ids.emplace_back(g_idx);
                values.emplace_back(_value);
            }
        }
    }

private:
    double _value;
    GeoLib::GeoObject const* const _geometry;
};


}   // namespace ProcessLib

#endif  // PROCESS_LIB_BOUNDARY_CONDITION_H_
