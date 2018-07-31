/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#ifdef OGS_USE_PYTHON

#include <pybind11/pybind11.h>
#include "BoundaryCondition.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/IndexValueVector.h"

#include "GenericNaturalBoundaryConditionLocalAssembler.h"

namespace ProcessLib
{
class PythonBoundaryConditionPythonSideInterface;

struct PythonBoundaryConditionData
{
    PythonBoundaryConditionPythonSideInterface* bc_object;
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk;
    std::size_t const bulk_mesh_id;
    int const global_component_id;

    const MeshLib::Mesh& boundary_mesh;
};

class PythonBoundaryCondition : public BoundaryCondition
{
public:
    PythonBoundaryCondition(PythonBoundaryConditionData&& bc_data,
                            unsigned const integration_order,
                            unsigned const shapefunction_order,
                            unsigned const global_dim,
                            bool const flush_stdout);

    void getEssentialBCValues(
        const double t, const GlobalVector& x,
        NumLib::IndexValueVector<GlobalIndexType>& bc_values) const override;

    void applyNaturalBC(const double t, const GlobalVector& x, GlobalMatrix& K,
                        GlobalVector& b, GlobalMatrix* Jac) override;

private:
    PythonBoundaryConditionData _bc_data;

    /// Local dof table, a subset of the global one restricted to the
    /// participating #_elements of the boundary condition.
    std::unique_ptr<NumLib::LocalToGlobalIndexMap> _dof_table_boundary;

    std::unique_ptr<MeshLib::MeshSubset const> _mesh_subset_all_nodes;

    /// Integration order for integration over the lower-dimensional elements
    unsigned const _integration_order;

    /// Local assemblers for each element of #_elements.
    std::vector<
        std::unique_ptr<GenericNaturalBoundaryConditionLocalAssemblerInterface>>
        _local_assemblers;

    bool const _flush_stdout;
};

std::unique_ptr<PythonBoundaryCondition> createPythonBoundaryCondition(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table, std::size_t bulk_mesh_id,
    int const variable_id, int const component_id, bool is_axially_symmetric,
    unsigned const integration_order, unsigned const shapefunction_order,
    unsigned const global_dim);

}  // namespace ProcessLib

#endif
