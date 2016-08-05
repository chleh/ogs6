/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "AbstractJacobianAssembler.h"

namespace ProcessLib
{
class CentralDifferencesJacobianAssembler : public AbstractJacobianAssembler
{
public:
    void assembleWithJacobian(
        LocalAssemblerInterface& local_assembler, double const t,
        std::vector<double> const& local_x,
        std::vector<double> const& local_xdot, const double dxdot_dx,
        const double dx_dx, std::vector<double>& local_M_data,
        std::vector<double>& local_K_data, std::vector<double>& local_b_data,
        std::vector<double>& local_Jac_data) override;

private:
    std::vector<double> _local_M_data;
    std::vector<double> _local_K_data;
    std::vector<double> _local_b_data;
};

}  // ProcessLib
