/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <functional>
#include <vector>

#include <pybind11/pybind11.h>

#include "ProcessLib/BoundaryCondition/ActiveBoundaryConditionCriterion.h"

namespace ProcessLib
{
class ActiveBCCriterionPython final : public ActiveBoundaryConditionCriterion
{
public:
    ActiveBCCriterionPython(pybind11::object&& scope,
                            std::string const& python_function)
        : _scope(std::move(scope)), _python_function(python_function)
    {
    }

    std::size_t getActiveBC(const std::size_t current_active_bc,
                            const GlobalVector& x, const double t,
                            const NumLib::LocalToGlobalIndexMap& dof_table,
                            const std::size_t mesh_id,
                            const std::vector<std::size_t>& node_ids) override;

private:
    pybind11::object _scope;
    std::string const _python_function;
};

std::unique_ptr<ActiveBCCriterionPython> createActiveBCCriterionTESPython(
    BaseLib::ConfigTree const& config);

}  // ProcessLib
