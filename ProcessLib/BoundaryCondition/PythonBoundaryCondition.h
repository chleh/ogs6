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
struct PythonBoundaryConditionData
{
    pybind11::object scope;
    std::string bc_object;
};

class PythonBoundaryCondition : public BoundaryCondition
{
public:
    PythonBoundaryCondition(PythonBoundaryConditionData&& bc_data,
                            std::vector<std::size_t>&& mesh_node_ids,
                            std::vector<MeshLib::Element*>&& elements,
                            const MeshLib::Mesh& mesh,
                            NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
                            int const variable_id,
                            int const component_id,
                            unsigned const integration_order,
                            unsigned const shapefunction_order);

    void getEssentialBCValues(
        const double t,
        NumLib::IndexValueVector<GlobalIndexType>& bc_values) const override;

    void applyNaturalBC(const double t,
                        const GlobalVector& x,
                        GlobalMatrix& K,
                        GlobalVector& b) override;

    ~PythonBoundaryCondition();

private:
    PythonBoundaryConditionData const _bc_data;

    std::vector<std::size_t> const _mesh_node_ids;
    std::vector<MeshLib::Element*> const _elements;
    const MeshLib::Mesh& _mesh;
    NumLib::LocalToGlobalIndexMap const& _dof_table_bulk;
    int const _variable_id;
    int const _component_id;

    std::unique_ptr<MeshLib::MeshSubset const> _mesh_subset_all_nodes;

    /// Local dof table, a subset of the global one restricted to the
    /// participating #_elements of the boundary condition.
    std::unique_ptr<NumLib::LocalToGlobalIndexMap> _dof_table_boundary;

    /// Integration order for integration over the lower-dimensional elements
    unsigned const _integration_order;

    /// Local assemblers for each element of #_elements.
    std::vector<
        std::unique_ptr<GenericNaturalBoundaryConditionLocalAssemblerInterface>>
        _local_assemblers;
};

std::unique_ptr<PythonBoundaryCondition> createPythonBoundaryCondition(
    BaseLib::ConfigTree const& config, std::vector<std::size_t>&& mesh_node_ids,
    std::vector<MeshLib::Element*>&& elements,
    NumLib::LocalToGlobalIndexMap const& dof_table, int const variable_id,
    int const component_id, MeshLib::Mesh const& mesh,
    unsigned const integration_order, unsigned const shapefunction_order,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters);

}  // namespace ProcessLib

#endif
