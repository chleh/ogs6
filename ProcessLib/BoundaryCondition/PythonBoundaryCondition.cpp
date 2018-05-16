/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "PythonBoundaryCondition.h"

#include <pybind11/eval.h>
#include <pybind11/stl.h>

#include "MeshLib/MeshSearch/NodeSearch.h"
#include "ProcessLib/Utils/CreateLocalAssemblers.h"
#include "ProcessLib/Utils/ProcessUtils.h"
#include "PythonBoundaryConditionLocalAssembler.h"

namespace ProcessLib
{
PythonBoundaryCondition::PythonBoundaryCondition(
    PythonBoundaryConditionData&& bc_data,
    std::vector<std::size_t>&& mesh_node_ids,
    std::vector<MeshLib::Element*>&& elements,
    int const variable_id,
    int const component_id,
    unsigned const integration_order,
    unsigned const shapefunction_order)
    : _bc_data(std::move(bc_data)),
      _mesh_node_ids(std::move(mesh_node_ids)),
      _elements(std::move(elements)),
      _variable_id(variable_id),
      _component_id(component_id),
      _integration_order(integration_order)
{
    if (_mesh_node_ids.empty())
        OGS_FATAL("no mesh node ids found");

    // check basic data consistency
    if (variable_id >=
            static_cast<int>(_bc_data.dof_table_bulk.getNumberOfVariables()) ||
        component_id >=
            _bc_data.dof_table_bulk.getNumberOfVariableComponents(variable_id))
    {
        OGS_FATAL(
            "Variable id or component id too high. Actual values: (%d, %d), "
            "maximum values: (%d, %d).",
            variable_id, component_id,
            _bc_data.dof_table_bulk.getNumberOfVariables(),
            _bc_data.dof_table_bulk.getNumberOfVariableComponents(variable_id));
    }

    std::vector<MeshLib::Node*> nodes = MeshLib::getUniqueNodes(_elements);
    DBUG("Found %d nodes for Natural BCs for the variable %d and component %d",
         nodes.size(), variable_id, component_id);

    auto const& mesh_subsets =
        _bc_data.dof_table_bulk.getMeshSubsets(variable_id, component_id);

    // TODO extend the node intersection to all parts of mesh_subsets, i.e.
    // to each of the MeshSubset in the mesh_subsets.
    _mesh_subset_all_nodes.reset(
        mesh_subsets.getMeshSubset(0).getIntersectionByNodes(nodes));
    MeshLib::MeshSubsets all_mesh_subsets{_mesh_subset_all_nodes.get()};

    // Create local DOF table from intersected mesh subsets for the given
    // variable and component ids.
    _bc_data.dof_table_boundary.reset(
        _bc_data.dof_table_bulk.deriveBoundaryConstrainedMap(
            variable_id, {component_id}, std::move(all_mesh_subsets),
            _elements));

    createLocalAssemblers<PythonConditionLocalAssembler>(
        _bc_data.mesh.getDimension(), _elements, *_bc_data.dof_table_boundary,
        shapefunction_order, _local_assemblers,
        _bc_data.mesh.isAxiallySymmetric(), _integration_order, _bc_data);
}

void PythonBoundaryCondition::getEssentialBCValues(
    const double t, NumLib::IndexValueVector<GlobalIndexType>& bc_values) const
{
    auto* bc = _bc_data.scope[_bc_data.bc_object.c_str()]
                   .cast<::PyBoundaryCondition*>();

    auto const mesh_id = _bc_data.mesh.getID();
    auto const nodes = _bc_data.mesh.getNodes();

    bc_values.ids.clear();
    bc_values.values.clear();

    bc_values.ids.reserve(_mesh_node_ids.size());
    bc_values.values.reserve(_mesh_node_ids.size());

    SpatialPosition pos;

    for (auto const id : _mesh_node_ids)
    {
        auto* xs = nodes[id]->getCoords();  // TODO DDC problems?
        auto val = bc->getDirichletBCValue(t, {xs[0], xs[1], xs[2]});
        if (!bc->isOverriddenEssential())
        {
            DBUG(
                "Method `getDirichletBCValue' not overridden in Python "
                "script.");
            return;
        }
        if (val.first)
        {
            pos.setNodeID(id);
            MeshLib::Location l(mesh_id, MeshLib::MeshItemType::Node, id);
            const auto g_idx = _bc_data.dof_table_bulk.getGlobalIndex(
                l, _variable_id, _component_id);
            if (g_idx == NumLib::MeshComponentMap::nop)
                continue;
            // For the DDC approach (e.g. with PETSc option), the negative
            // index of g_idx means that the entry by that index is a ghost
            // one, which should be dropped. Especially for PETSc routines
            // MatZeroRows and MatZeroRowsColumns, which are called to apply
            // the Dirichlet BC, the negative index is not accepted like
            // other matrix or vector PETSc routines. Therefore, the
            // following if-condition is applied.
            if (g_idx >= 0)
            {
                bc_values.ids.emplace_back(g_idx);
                bc_values.values.emplace_back(val.second);
            }
        }
    }
}

void PythonBoundaryCondition::applyNaturalBC(const double t,
                                             const GlobalVector& x,
                                             GlobalMatrix& K,
                                             GlobalVector& b)
{
    try
    {
        GlobalExecutor::executeMemberOnDereferenced(
            &GenericNaturalBoundaryConditionLocalAssemblerInterface::assemble,
            _local_assemblers, *_bc_data.dof_table_boundary, t, x, K, b);
    }
    catch (::PyNotOverridden const& /*e*/)
    {
        DBUG("Method `getFlux' not overridden in Python script.");
    }
}

PythonBoundaryCondition::~PythonBoundaryCondition()
{
    for (auto e : _elements)
        delete e;
}

std::unique_ptr<PythonBoundaryCondition> createPythonBoundaryCondition(
    BaseLib::ConfigTree const& config, std::vector<std::size_t>&& mesh_node_ids,
    std::vector<MeshLib::Element*>&& elements,
    NumLib::LocalToGlobalIndexMap const& dof_table, int const variable_id,
    int const component_id, const MeshLib::Mesh& mesh,
    unsigned const integration_order, unsigned const shapefunction_order,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters)
{
    config.checkConfigParameter("type", "Python");

    // //// Parse Python //////////////////////

    auto const script = config.getConfigParameter<std::string>("python_script");
    auto const bc_object = config.getConfigParameter<std::string>("bc_object");

    namespace py = pybind11;
    // using namespace py::literals;

    // Evaluate in scope of main module
    py::module module = py::module::import("__main__");
    py::object scope = module.attr("__dict__");

    // TODO
    // http://pybind11.readthedocs.io/en/stable/advanced/embedding.html#adding-embedded-modules
    // first attempt on that failed
    py::module ogs = module.def_submodule("OpenGeoSys", "NO HELP AVAILABLE");

    py::class_<::PyBoundaryCondition, ::PyBoundaryConditionImpl> pybc(
        ogs, "BoundaryCondition");
    pybc.def(py::init());
    pybc.def("getDirichletBCValue",
             &::PyBoundaryCondition::getDirichletBCValue);
    pybc.def("getFlux", &::PyBoundaryCondition::getFlux);

    py::eval_file(script, scope);

    if (!scope.contains(bc_object))
        OGS_FATAL("Function `%s' is not defined in file `%s'.",
                  bc_object.c_str(), script.c_str());

    return std::make_unique<PythonBoundaryCondition>(
        PythonBoundaryConditionData{std::move(scope), bc_object, dof_table,
                                    nullptr, mesh},
        std::move(mesh_node_ids), std::move(elements), variable_id,
        component_id, integration_order, shapefunction_order);
}

}  // namespace ProcessLib
