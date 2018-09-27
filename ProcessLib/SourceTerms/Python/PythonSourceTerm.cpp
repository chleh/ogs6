/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "PythonSourceTerm.h"

#include <pybind11/pybind11.h>
#include <iostream>

#include "MeshLib/MeshSearch/NodeSearch.h"
#include "ProcessLib/Utils/CreateLocalAssemblers.h"
#include "ProcessLib/Utils/ProcessUtils.h"
#include "PythonSourceTermLocalAssembler.h"

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
    //! To flush or not to flush.
    const bool _flush;
};
}  // anonymous namespace

namespace ProcessLib
{
PythonSourceTerm::PythonSourceTerm(
    NumLib::LocalToGlobalIndexMap const& source_term_dof_table,
    PythonSourceTermData&& source_term_data, unsigned const integration_order,
    unsigned const shapefunction_order, unsigned const global_dim,
    bool const flush_stdout)
    : SourceTerm(source_term_dof_table),
      _source_term_data(std::move(source_term_data)),
      _flush_stdout(flush_stdout)
{
    std::vector<MeshLib::Node*> const& source_term_nodes =
        _source_term_data.source_term_mesh.getNodes();
    MeshLib::MeshSubset source_term_mesh_subset(
        _source_term_data.source_term_mesh, source_term_nodes);

    // Create local DOF table from the source term mesh subset for the given
    // variable and component id.
    _dof_table_source_term =
        _source_term_data.dof_table_bulk.deriveBoundaryConstrainedMap(
            std::move(source_term_mesh_subset));

    createLocalAssemblers<PythonSourceTermLocalAssembler>(
        global_dim, _source_term_data.source_term_mesh.getElements(),
        *_dof_table_source_term, shapefunction_order, _local_assemblers,
        _source_term_data.source_term_mesh.isAxiallySymmetric(),
        integration_order, _source_term_data);
}

void PythonSourceTerm::integrate(const double t, const GlobalVector& x,
                                 GlobalVector& b, GlobalMatrix* Jac) const
{
    FlushStdoutGuard guard(_flush_stdout);

    try
    {
        GlobalExecutor::executeMemberOnDereferenced(
            &PythonSourceTermLocalAssemblerInterface::assemble,
            _local_assemblers, *_dof_table_source_term, t, x, b, Jac);
    }
    catch (MethodNotOverriddenInDerivedClassException const& /*e*/)
    {
        DBUG("Method `getFlux' not overridden in Python script.");
    }
}

std::unique_ptr<PythonSourceTerm> createPythonSourceTerm(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& source_term_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table, std::size_t bulk_mesh_id,
    int const variable_id, int const component_id,
    unsigned const integration_order, unsigned const shapefunction_order,
    unsigned const global_dim)
{
    //! \ogs_file_param{prj__process_variables__process_variable__source_term__source_term__type}
    config.checkConfigParameter("type", "Python");

    //! \ogs_file_param{prj__process_variables__process_variable__source_term__source_term__Python__source_term_object}
    auto const source_term_object = config.getConfigParameter<std::string>("source_term_object");
    //! \ogs_file_param{prj__process_variables__process_variable__source_term__source_term__Python__flush_stdout}
    auto const flush_stdout = config.getConfigParameter("flush_stdout", false);

    // Evaluate Python code in scope of main module
    pybind11::object scope =
        pybind11::module::import("__main__").attr("__dict__");

    if (!scope.contains(source_term_object))
        OGS_FATAL(
            "Function `%s' is not defined in the python script file, or there "
            "was no python script file specified.",
            source_term_object.c_str());

    auto* source_term = scope[source_term_object.c_str()]
                   .cast<PythonSourceTermPythonSideInterface*>();

    if (variable_id >= static_cast<int>(dof_table.getNumberOfVariables()) ||
        component_id >= dof_table.getNumberOfVariableComponents(variable_id))
    {
        OGS_FATAL(
            "Variable id or component id too high. Actual values: (%d, %d), "
            "maximum values: (%d, %d).",
            variable_id, component_id, dof_table.getNumberOfVariables(),
            dof_table.getNumberOfVariableComponents(variable_id));
    }

    // In case of partitioned mesh the source_term could be empty, i.e. there is no
    // source_term condition.
#ifdef USE_PETSC
    // This can be extracted to createSourceTerm() but then the config
    // parameters are not read and will cause an error.
    // TODO (naumov): Add a function to ConfigTree for skipping the tags of the
    // subtree and move the code up in createSourceTerm().
    if (source_term_mesh.getDimension() == 0 &&
        source_term_mesh.getNumberOfNodes() == 0 &&
        source_term_mesh.getNumberOfElements() == 0)
    {
        return nullptr;
    }
#endif  // USE_PETSC

    return std::make_unique<PythonSourceTerm>(
        PythonSourceTermData{
            source_term, dof_table, bulk_mesh_id,
            dof_table.getGlobalComponent(variable_id, component_id),
            source_term_mesh},
        integration_order, shapefunction_order, global_dim, flush_stdout);
}

}  // namespace ProcessLib
