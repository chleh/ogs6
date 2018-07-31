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
//! Optionally flushes the standard output upon creation and destruction.
//! Can be used to improve the debug output readability when printing debug
//! messages both from OGS and from Python.
class FlushStdoutGuard
{
public:
    //! Optionally flushes C++ stdout before running Python code.
    explicit FlushStdoutGuard(bool const flush) : _flush(flush)
    {
        if (!flush)
            return;

        LOGOG_COUT << std::flush;
    }

    //! Optionally flushes Python's stdout after running Python code.
    ~FlushStdoutGuard()
    {
        if (!_flush)
            return;

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
    unsigned const integration_order,
    unsigned const shapefunction_order,
    unsigned const global_dim,
    bool const flush_stdout)
    : _bc_data(std::move(bc_data)),
      _integration_order(integration_order),
      _flush_stdout(flush_stdout)
{
    std::vector<MeshLib::Node*> const& bc_nodes =
        _bc_data.boundary_mesh.getNodes();
    MeshLib::MeshSubset bc_mesh_subset(_bc_data.boundary_mesh, bc_nodes);

    // Create local DOF table from the bc mesh subset for the given variable and
    // component id.
    _dof_table_boundary = _bc_data.dof_table_bulk.deriveBoundaryConstrainedMap(
        std::move(bc_mesh_subset));

    createLocalAssemblers<PythonBoundaryConditionLocalAssembler>(
        global_dim, _bc_data.boundary_mesh.getElements(), *_dof_table_boundary,
        shapefunction_order, _local_assemblers,
        _bc_data.boundary_mesh.isAxiallySymmetric(), _integration_order,
        _bc_data);
}

void PythonBoundaryCondition::getEssentialBCValues(
    const double t, GlobalVector const& x,
    NumLib::IndexValueVector<GlobalIndexType>& bc_values) const
{
    FlushStdoutGuard guard(_flush_stdout);

    auto const nodes = _bc_data.boundary_mesh.getNodes();

    auto const& bulk_node_ids_map =
        *_bc_data.boundary_mesh.getProperties().getPropertyVector<std::size_t>(
            "bulk_node_ids");

    bc_values.ids.clear();
    bc_values.values.clear();

    bc_values.ids.reserve(_bc_data.boundary_mesh.getNumberOfNodes());
    bc_values.values.reserve(_bc_data.boundary_mesh.getNumberOfNodes());

    std::vector<double> primary_variables;

    for (auto const* node : _bc_data.boundary_mesh.getNodes())
    {
        auto const boundary_node_id = node->getID();
        auto const bulk_node_id = bulk_node_ids_map[boundary_node_id];

        // gather primary variable values
        primary_variables.clear();
        auto const num_var = _dof_table_boundary->getNumberOfVariables();
        for (int var = 0; var < num_var; ++var)
        {
            auto const num_comp =
                _dof_table_boundary->getNumberOfVariableComponents(var);
            for (int comp = 0; comp < num_comp; ++comp)
            {
                MeshLib::Location loc{_bc_data.bulk_mesh_id,
                                      MeshLib::MeshItemType::Node,
                                      bulk_node_id};
                auto const dof_idx =
                    _bc_data.dof_table_bulk.getGlobalIndex(loc, var, comp);
                primary_variables.push_back(x[dof_idx]);
            }
        }

        auto* xs = nodes[boundary_node_id]->getCoords();  // TODO DDC problems?
        auto val = _bc_data.bc_object->getDirichletBCValue(
            t, {xs[0], xs[1], xs[2]}, boundary_node_id, primary_variables);
        if (!_bc_data.bc_object->isOverriddenEssential())
        {
            DBUG(
                "Method `getDirichletBCValue' not overridden in Python "
                "script.");
            return;
        }
        if (val.first)
        {
            MeshLib::Location l(_bc_data.bulk_mesh_id,
                                MeshLib::MeshItemType::Node, bulk_node_id);
            const auto dof_idx = _bc_data.dof_table_bulk.getGlobalIndex(
                l, _bc_data.global_component_id);
            if (dof_idx == NumLib::MeshComponentMap::nop)
                continue;
            // For the DDC approach (e.g. with PETSc option), the negative
            // index of g_idx means that the entry by that index is a ghost
            // one, which should be dropped. Especially for PETSc routines
            // MatZeroRows and MatZeroRowsColumns, which are called to apply
            // the Dirichlet BC, the negative index is not accepted like
            // other matrix or vector PETSc routines. Therefore, the
            // following if-condition is applied.
            if (dof_idx >= 0)
            {
                bc_values.ids.emplace_back(dof_idx);
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
    catch (MethodNotOverriddenInDerivedClassException const& /*e*/)
    {
        DBUG("Method `getFlux' not overridden in Python script.");
    }
}

std::unique_ptr<PythonBoundaryCondition> createPythonBoundaryCondition(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table, std::size_t bulk_mesh_id,
    int const variable_id, int const component_id, bool is_axially_symmetric,
    unsigned const integration_order, unsigned const shapefunction_order,
    unsigned const global_dim)
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

    auto* bc = scope[bc_object.c_str()]
                   .cast<PythonBoundaryConditionPythonSideInterface*>();

    if (variable_id >= static_cast<int>(dof_table.getNumberOfVariables()) ||
        component_id >= dof_table.getNumberOfVariableComponents(variable_id))
    {
        OGS_FATAL(
            "Variable id or component id too high. Actual values: (%d, %d), "
            "maximum values: (%d, %d).",
            variable_id, component_id, dof_table.getNumberOfVariables(),
            dof_table.getNumberOfVariableComponents(variable_id));
    }

    return std::make_unique<PythonBoundaryCondition>(
        PythonBoundaryConditionData{
            bc, dof_table, bulk_mesh_id,
            dof_table.getGlobalComponent(variable_id, component_id), bc_mesh},
        integration_order, shapefunction_order, global_dim, flush_stdout);
}

}  // namespace ProcessLib
