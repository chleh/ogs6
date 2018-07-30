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
#include <iostream>

#include "MeshLib/MeshSearch/NodeSearch.h"
#include "ProcessLib/Utils/CreateLocalAssemblers.h"
#include "ProcessLib/Utils/ProcessUtils.h"
#include "PythonBoundaryConditionLocalAssembler.h"

namespace
{
class FlushStdoutGuard
{
public:
    explicit FlushStdoutGuard(bool const flush) : _flush(flush)
    {
        // flush std::cout before running Python code
        if (!flush)
            return;

        // std::cout << std::flush;
        LOGOG_COUT << std::flush;
    }

    ~FlushStdoutGuard()
    {
        if (!_flush)
            return;

        // flush Python's stdout after running python code
        using namespace pybind11::literals;
        pybind11::print("end"_a = "", "flush"_a = true);
    }

private:
    const bool _flush;
};
}  // anonymous namespace

namespace ProcessLib
{
PythonBoundaryCondition::PythonBoundaryCondition(
    PythonBoundaryConditionData&& bc_data,
    std::vector<std::size_t>&& mesh_node_ids,
    std::vector<MeshLib::Element*>&& elements,
    unsigned const integration_order,
    unsigned const shapefunction_order,
    bool const flush_stdout)
    : _bc_data(std::move(bc_data)),
      _mesh_node_ids(std::move(mesh_node_ids)),
      _elements(std::move(elements)),
      _integration_order(integration_order),
      _flush_stdout(flush_stdout)
{
    if (_mesh_node_ids.empty())
        OGS_FATAL("no mesh node ids found");

    std::vector<MeshLib::Node*> nodes = MeshLib::getUniqueNodes(_elements);

    // FIXME
#if 0
    auto const& mesh_subset =
        _bc_data.dof_table_bulk.getMeshSubset(_bc_data.global_component_id);

    mesh_subset.getI

    // TODO extend the node intersection to all parts of mesh_subsets, i.e.
    // to each of the MeshSubset in the mesh_subsets.
    _mesh_subset_all_nodes.reset(
        mesh_subsets.getMeshSubset(0).getIntersectionByNodes(nodes));
    MeshLib::MeshSubsets all_mesh_subsets{_mesh_subset_all_nodes.get()};

    // Create local DOF table from intersected mesh subsets for the given
    // variable and component ids.
    _dof_table_boundary = _bc_data.dof_table_bulk.deriveBoundaryConstrainedMap(
        std::move(all_mesh_subsets), _elements);

    createLocalAssemblers<PythonBoundaryConditionLocalAssembler>(
        _bc_data.mesh.getDimension(), _elements, *_dof_table_boundary,
        shapefunction_order, _local_assemblers,
        _bc_data.mesh.isAxiallySymmetric(), _integration_order, _bc_data);
#endif
}

void PythonBoundaryCondition::getEssentialBCValues(
    const double t, GlobalVector const& x,
    NumLib::IndexValueVector<GlobalIndexType>& bc_values) const
{
    FlushStdoutGuard guard(_flush_stdout);

    auto const mesh_id = _bc_data.mesh.getID();
    auto const nodes = _bc_data.mesh.getNodes();

    bc_values.ids.clear();
    bc_values.values.clear();

    bc_values.ids.reserve(_mesh_node_ids.size());
    bc_values.values.reserve(_mesh_node_ids.size());

    SpatialPosition pos;
    std::vector<double> primary_variables;

    for (auto const node_id : _mesh_node_ids)
    {
        // gather primary variable values
        primary_variables.clear();
        auto const num_var = _bc_data.dof_table_bulk.getNumberOfVariables();
        for (int var = 0; var < num_var; ++var)
        {
            auto const num_comp =
                _bc_data.dof_table_bulk.getNumberOfVariableComponents(var);
            for (int comp = 0; comp < num_comp; ++comp)
            {
                MeshLib::Location loc{_bc_data.mesh.getID(),
                                      MeshLib::MeshItemType::Node, node_id};
                auto const dof_idx =
                    _bc_data.dof_table_bulk.getGlobalIndex(loc, var, comp);
                primary_variables.push_back(x[dof_idx]);
            }
        }

        auto* xs = nodes[node_id]->getCoords();  // TODO DDC problems?
        auto val = _bc_data.bc_object->getDirichletBCValue(
            t, {xs[0], xs[1], xs[2]}, node_id, primary_variables);
        if (!_bc_data.bc_object->isOverriddenEssential())
        {
            DBUG(
                "Method `getDirichletBCValue' not overridden in Python "
                "script.");
            return;
        }
        if (val.first)
        {
            pos.setNodeID(node_id);
            MeshLib::Location l(mesh_id, MeshLib::MeshItemType::Node, node_id);
            const auto g_idx = _bc_data.dof_table_bulk.getGlobalIndex(
                l, _bc_data.global_component_id);
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
                                             GlobalMatrix& K, GlobalVector& b,
                                             GlobalMatrix* Jac)
{
    FlushStdoutGuard guard(_flush_stdout);

    try
    {
        GlobalExecutor::executeMemberOnDereferenced(
            &GenericNaturalBoundaryConditionLocalAssemblerInterface::assemble,
            _local_assemblers, *_dof_table_boundary, t, x, K, b, Jac);
    }
    catch (PyNotOverridden const& /*e*/)
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
    unsigned const integration_order, unsigned const shapefunction_order)
{
    config.checkConfigParameter("type", "Python");

    auto const bc_object = config.getConfigParameter<std::string>("bc_object");
    auto const flush_stdout = config.getConfigParameter("flush_stdout", false);

    // Evaluate in scope of main module
    pybind11::object scope =
        pybind11::module::import("__main__").attr("__dict__");

    if (!scope.contains(bc_object))
        OGS_FATAL(
            "Function `%s' is not defined in the python script file, or there "
            "was no python script file specified.",
            bc_object.c_str());

    auto* bc = scope[bc_object.c_str()].cast<PyBoundaryCondition*>();

    return std::make_unique<PythonBoundaryCondition>(
        PythonBoundaryConditionData{
            bc, dof_table,
            dof_table.getGlobalComponent(variable_id, component_id), mesh},
        std::move(mesh_node_ids), std::move(elements), integration_order,
        shapefunction_order, flush_stdout);
}

}  // namespace ProcessLib
