/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CentralDifferencesJacobianAssembler.h"

namespace ProcessLib
{
void CentralDifferencesJacobianAssembler::assembleWithJacobian(
    LocalAssemblerInterface& local_assembler, const double t,
    const std::vector<double>& local_x, const std::vector<double>& local_xdot,
    const double dxdot_dx, const double dx_dx,
    std::vector<double>& local_M_data, std::vector<double>& local_K_data,
    std::vector<double>& local_b_data,
    std::vector<double>& local_Jac_data)
{
    // TODO implement
}

}  // ProcessLib
