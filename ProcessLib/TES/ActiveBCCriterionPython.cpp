/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifdef OGS_USE_PYTHON

#include "ActiveBCCriterionPython.h"
#include <pybind11/eval.h>
#include <pybind11/stl.h>
#include "NumLib/DOF/LocalToGlobalIndexMap.h"

namespace ProcessLib
{
std::size_t ActiveBCCriterionPython::getActiveBC(
    const std::size_t current_active_bc, const GlobalVector& x, const double t,
    const NumLib::LocalToGlobalIndexMap& dof_table, const std::size_t mesh_id,
    const std::vector<std::size_t>& node_ids)
{
    if (node_ids.empty())
        OGS_FATAL("No node for measuring provided.");

    auto const num_dofs = 3;
    if (dof_table.getNumberOfComponents() != num_dofs)
        OGS_FATAL("Wrong number of DOFs.");

    MeshLib::Location loc(mesh_id, MeshLib::MeshItemType::Node, 0);

    auto const num_nodes = node_ids.size();
    std::array<std::vector<double>, num_dofs> values;
    for (auto& v : values)
        v.resize(num_nodes);

    for (std::size_t n = 0; n < num_nodes; ++n)
    {
        loc.item_id = node_ids[n];
        for (int c = 0; c < num_dofs; ++c)
        {
            auto const dof_idx = dof_table.getLocalIndex(
                loc, c, x.getRangeBegin(), x.getRangeEnd());
            values[c][n] = x.get(dof_idx);
        }
    }

    return _scope[_python_function.c_str()](current_active_bc, t, values[0],
                                            values[1], values[2])
        .cast<std::size_t>();
}

std::unique_ptr<ActiveBCCriterionPython> createActiveBCCriterionTESPython(
    BaseLib::ConfigTree const& config)
{
    config.checkConfigParameter("type", "TESPython");

    auto const script = config.getConfigParameter<std::string>("python_script");
    auto const function =
        config.getConfigParameter<std::string>("python_function");

    namespace py = pybind11;

    // Evaluate in scope of main module
    py::module module = py::module::import("__main__");
    py::object scope = module.attr("__dict__");

    py::eval_file(script, scope);

    if (!scope.contains(function))
        OGS_FATAL("Function `%s' is not defined in file `%s'.",
                  function.c_str(), script.c_str());

    return std::unique_ptr<ActiveBCCriterionPython>(
        new ActiveBCCriterionPython(std::move(scope), function));
}

}  // ProcessLib

#endif
