/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <vector>

namespace MeshLib
{
class FEMMesh;
}

namespace NumLib
{
class AbstractDOFTable;
}  // namespace NumLib

namespace ProcessLib
{
class BoundaryCondition;
struct BoundaryConditionConfig;
struct ParameterBase;
class Process;

std::unique_ptr<BoundaryCondition> createBoundaryCondition(
    const BoundaryConditionConfig& config,
    const NumLib::AbstractDOFTable& dof_table,
    const MeshLib::FEMMesh& bulk_mesh, const int variable_id,
    const unsigned integration_order, const unsigned shapefunction_order,
    const std::vector<std::unique_ptr<ProcessLib::ParameterBase>>& parameters,
    const Process& process);

}  // namespace ProcessLib
