/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "UniformDirichletBoundaryCondition.h"

#include <algorithm>
#include <vector>

#include <logog/include/logog.hpp>

namespace ProcessLib
{
UniformDirichletBoundaryCondition::UniformDirichletBoundaryCondition(
    GeoLib::GeoObject const& geometry, BaseLib::ConfigTree const& config)
    : _geometry(geometry)
{
    DBUG("Constructing UniformDirichletBoundaryCondition from config.");
    //! \ogs_file_param{boundary_condition__type}
    config.checkConfigParameter("type", "UniformDirichlet");

    //! \ogs_file_param{boundary_condition__UniformDirichlet__value}
    _value = config.getConfigParameter<double>("value");
    DBUG("Using value %g", _value);
}

UniformDirichletBoundaryCondition::UniformDirichletBoundaryCondition(
    GeoLib::GeoObject const& geometry, double value)
    : _value(value), _geometry(geometry)
{
    DBUG("Constructed UniformDirichletBoundaryCondition using value %g",
         _value);
}

std::unique_ptr<DirichletBoundaryCondition>
UniformDirichletBoundaryCondition::getDirichletBoundaryCondition(
    MeshGeoToolsLib::MeshNodeSearcher& searcher,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    int const variable_id,
    int const component_id)
{
    // Find nodes' ids on the given mesh on which this boundary condition
    // is defined.
    std::vector<std::size_t> ids = searcher.getMeshNodeIDs(_geometry);

    NumLib::IndexValueVector<GlobalIndexType> bc;

    // convert mesh node ids to global index for the given component
    bc.ids.reserve(bc.ids.size() + ids.size());
    bc.values.reserve(bc.values.size() + ids.size());
    for (auto& id : ids)
    {
        MeshLib::Location l(searcher.getMeshId(), MeshLib::MeshItemType::Node,
                            id);
        // TODO: that might be slow, but only done once
        const auto g_idx =
            dof_table.getGlobalIndex(l, variable_id, component_id);
        // For the DDC approach (e.g. with PETSc option), the negative
        // index of g_idx means that the entry by that index is a ghost one,
        // which should be dropped. Especially for PETSc routines MatZeroRows
        // and MatZeroRowsColumns, which are called to apply the Dirichlet BC,
        // the negative index is not accepted like other matrix or vector
        // PETSc routines. Therefore, the following if-condition is applied.
        if (g_idx >= 0)
        {
            bc.ids.emplace_back(g_idx);
            bc.values.emplace_back(_value);
        }
    }

    return std::unique_ptr<DirichletBoundaryCondition>(
        new DirichletBoundaryCondition(std::move(bc)));
}

}   // namespace ProcessLib
