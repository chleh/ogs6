/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <map>
#include <unordered_set>

#include "MeshLib/DUNEMesh.h"

#include "BoundaryCondition.h"

#include "MeshLib/PropertyVector.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/IndexValueVector.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Utils/ProcessUtils.h"

namespace ProcessLib
{
template <int GlobalDim>
class NonuniformDirichletBoundaryConditionDUNE final : public BoundaryCondition
{
public:
    NonuniformDirichletBoundaryConditionDUNE(
        MeshLib::Mesh const& boundary_mesh,
        MeshLib::PropertyVector<double> const& values,
        MeshLib::DUNEMesh<GlobalDim> const& bulk_mesh,
        int const variable_id_bulk,
        int const component_id_bulk);

    void getEssentialBCValues(
        const double /*t*/, GlobalVector const& /*x*/,
        NumLib::IndexValueVector<GlobalIndexType>& bc_values) const override;

private:
    void onPostRefine(MeshLib::DUNEMesh<GlobalDim> const& mesh,
                      bool globally_refined,
                      const MeshLib::DUNEIdToIdxMappings&);

    //! Indices of vertices where Dirichlet BCs are defined.
    //!
    //! This is merely a cache for _map_boundary_vertex_id_to_value.
    std::vector<typename BaseLib::DUNEGridType<
        GlobalDim>::LeafGridView::IndexSet::IndexType>
        _dirichlet_vertex_indices;

    //! Dirichlet BC Values corresponding to the specified vertex indices.
    //!
    //! This is merely a cache for _map_boundary_vertex_id_to_value.
    std::vector<double> _dirichlet_values;

    // int const _variable_id_bulk; // TODO [DUNE] implement
    int const _component_id_bulk;  //! Component to which this BC applies.

    //! Maps DUNE's persistent ids of the Dirichlet vertices to the Dirichlet
    //! values.
    std::map<typename BaseLib::DUNEGridType<GlobalDim>::LocalIdSet::IdType,
             double>
        _map_boundary_vertex_id_to_value;

    // TODO [DUNE] make const
    //! The set of all DUNE boundary segment indices belonging to this BC.
    std::unordered_set<std::size_t> _boundary_segment_indices;
};

std::unique_ptr<ProcessLib::BoundaryCondition>
createNonuniformDirichletBoundaryCondition(
    const ProcessLib::BoundaryConditionConfig& config,
    const NumLib::AbstractDOFTable& /*dof_table*/,
    const MeshLib::FEMMesh& bulk_mesh, const int variable_id);
}  // namespace ProcessLib
