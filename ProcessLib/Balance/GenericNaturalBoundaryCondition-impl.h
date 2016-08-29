/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "GenericNaturalBoundaryCondition.h"
#include "ProcessLib/Utils/CreateLocalAssemblers.h"
#include "MeshLib/MeshSearch/NodeSearch.h"
#include "GenericNaturalBoundaryConditionLocalAssembler.h"

namespace ProcessLib
{
BalanceProcess::BalanceProcess(
    unsigned const integration_order,
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk, int const variable_id,
    int const component_id, unsigned const global_dim,
    std::vector<MeshLib::Element*>&& elements, ...)
    : // ...,
      _elements(std::move(elements)), // or surface mesh
      _integration_order(integration_order)
{
    assert(component_id <
           static_cast<int>(dof_table_bulk.getNumberOfComponents()));

    // boundary dof table not needed, because you write per-element values
#if 0
    std::vector<MeshLib::Node*> nodes = MeshLib::getUniqueNodes(_elements);

    auto const& mesh_subsets =
        dof_table_bulk.getMeshSubsets(variable_id, component_id);

    // TODO extend the node intersection to all parts of mesh_subsets, i.e.
    // to each of the MeshSubset in the mesh_subsets.
    _mesh_subset_all_nodes.reset(
        mesh_subsets.getMeshSubset(0).getIntersectionByNodes(nodes));
    std::unique_ptr<MeshLib::MeshSubsets> all_mesh_subsets{
        new MeshLib::MeshSubsets{_mesh_subset_all_nodes.get()}};

    // Create local DOF table from intersected mesh subsets for the given
    // variable and component ids.
    _dof_table_boundary.reset(dof_table_bulk.deriveBoundaryConstrainedMap(
        variable_id, component_id, std::move(all_mesh_subsets), _elements));
#endif

    createLocalAssemblers<LocalAssemblerImplementation>(
        global_dim, _elements, _integration_order, _local_assemblers,

                _elements, map_surface_to_bulk // or surface mesh
                );
}

template <typename BoundaryConditionData,
          template <typename, typename, unsigned>
          class LocalAssemblerImplementation>
BalanceProcess<
    BoundaryConditionData,
    LocalAssemblerImplementation>::~BalanceProcess()
{
    for (auto e : _elements)
        delete e;
}

template <typename BoundaryConditionData,
          template <typename, typename, unsigned>
          class LocalAssemblerImplementation>
void BalanceProcess<
    BoundaryConditionData,
    LocalAssemblerImplementation>::apply(const double t,
                                         const GlobalVector& x,
                                         GlobalMatrix& K,
                                         GlobalVector& b)
{
    GlobalExecutor::executeMemberOnDereferenced(
        &BalanceProcessLocalAssemblerInterface::assemble,
        _local_assemblers, *_dof_table_boundary, t, x, K, b);
}

}  // ProcessLib
