/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CentralDifferencesJacobianAssembler.h"
#include <cassert>
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "LocalAssemblerInterface.h"

namespace ProcessLib
{
void CentralDifferencesJacobianAssembler::assembleWithJacobian(
    LocalAssemblerInterface& local_assembler, const double t,
    const std::vector<double>& local_x_data,
    const std::vector<double>& local_xdot_data, const double dxdot_dx,
    const double dx_dx, std::vector<double>& local_M_data,
    std::vector<double>& local_K_data, std::vector<double>& local_b_data,
    std::vector<double>& local_Jac_data)
{
    auto const eps = 1e-8;
    auto const num_r_c = static_cast<Eigen::Index>(local_x_data.size());

    auto const local_x = MathLib::toVector(local_x_data, num_r_c);
    auto const local_xdot = MathLib::toVector(local_xdot_data, num_r_c);

    auto local_Jac = MathLib::toZeroedMatrix(local_Jac_data, num_r_c, num_r_c);
    _local_x_perturbed_data = local_x_data;

    // Residual  res := M xdot + K x - b
    // Computing Jac := dres/dx
    //                = M dxdot/dx + dM/dx xdot + K dx/dx + dK/dx x - db/dx
    //                  (Note: dM/dx and dK/dx actually have the second and
    //                  third index transposed.)
    // The loop computes the dM/dx, dK/dx and db/dx terms, the rest is computed
    // afterwards.
    for (Eigen::Index i = 0; i < num_r_c; ++i)
    {
        _local_x_perturbed_data[i] += eps;
        local_assembler.assemble(t, _local_x_perturbed_data, local_M_data,
                                 local_K_data, local_b_data);

        _local_x_perturbed_data[i] = local_x_data[i] - eps;
        local_assembler.assemble(t, _local_x_perturbed_data, _local_M_data,
                                 _local_K_data, _local_b_data);

        _local_x_perturbed_data[i] = local_x_data[i];

        if (!local_M_data.empty()) {
            auto const local_M_p =
                MathLib::toMatrix(local_M_data, num_r_c, num_r_c);
            auto const local_M_m =
                MathLib::toMatrix(_local_M_data, num_r_c, num_r_c);
            local_Jac.col(i).noalias() +=
                // dM/dxi * x_dot
                (local_M_p - local_M_m) * local_xdot / (2.0 * eps);
            local_M_data.clear();
            _local_M_data.clear();
        }
        if (!local_K_data.empty()) {
            auto const local_K_p =
                MathLib::toMatrix(local_K_data, num_r_c, num_r_c);
            auto const local_K_m =
                MathLib::toMatrix(_local_K_data, num_r_c, num_r_c);
            local_Jac.col(i).noalias() +=
                // dK/dxi * x
                (local_K_p - local_K_m) * local_x / (2.0 * eps);
            local_K_data.clear();
            _local_K_data.clear();
        }
        if (!local_b_data.empty()) {
            auto const local_b_p = MathLib::toVector(local_b_data, num_r_c);
            auto const local_b_m = MathLib::toVector(_local_b_data, num_r_c);
            local_Jac.col(i).noalias() -=
                // db/dxi
                (local_b_p - local_b_m) / (2.0 * eps);
            local_b_data.clear();
            _local_b_data.clear();
        }
    }

    // Assemble with unperturbed local x.
    local_assembler.assemble(t, local_x_data, local_M_data, local_K_data,
                             local_b_data);

    // Compute remaining terms of the Jacobian.
    if (dxdot_dx != 0 && !local_M_data.empty()) {
        auto local_M = MathLib::toMatrix(local_M_data, num_r_c, num_r_c);
        local_Jac.noalias() += local_M * dxdot_dx;
    }
    if (dx_dx != 0 && !local_K_data.empty()) {
        auto local_K = MathLib::toMatrix(local_K_data, num_r_c, num_r_c);
        local_Jac.noalias() += local_K * dx_dx;
    }
}

}  // ProcessLib
