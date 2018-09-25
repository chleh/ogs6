/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>
#include <memory>

#include "ProcessLib/Parameter/Parameter.h"

namespace MeshLib
{
class Mesh;
}

namespace NumLib
{
class AbstractDOFTable;
}  // namespace NumLib

namespace ProcessLib
{
class SourceTerm;
struct SourceTermConfig;
}  // namespace ProcessLib

namespace ProcessLib
{
std::unique_ptr<SourceTerm> createSourceTerm(
    const SourceTermConfig& config, const NumLib::AbstractDOFTable& dof_table,
    const MeshLib::FEMMesh& mesh, const int variable_id,
    const unsigned integration_order, const unsigned shapefunction_order,
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters);

}  // namespace ProcessLib
