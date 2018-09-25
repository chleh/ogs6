/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "BaseLib/DUNEConfig.h"

#if OGS_USE_DUNE

#include "ProcessLib/ProcessVariable.h"
#include "ProcessOutput.h"
#include "SecondaryVariable.h"

namespace ProcessLib
{
//! Writes output to the given \c file_name using the VTU file format.
///
/// See Output::_output_file_data_mode documentation for the data_mode
/// parameter.
void doProcessOutputDUNE(
    std::string const& file_name,
    bool const compress_output,
    int const data_mode,
    const double t,
    GlobalVector const& x,
    MeshLib::FEMMesh& fem_mesh,
    NumLib::AbstractDOFTable const& dof_table,
    std::vector<std::reference_wrapper<ProcessVariable>> const&
        process_variables,
    SecondaryVariableCollection secondary_variables,
    ProcessOutput const& process_output);

}  // namespace ProcessLib

#endif
