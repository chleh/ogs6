/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CentralDifferencesJacobianAssembler.h"
#include <Eigen/Core>
#include <cassert>
#include "LocalAssemblerInterface.h"

Eigen::Map<
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
toMatrix(std::vector<double>& data, std::size_t num_r_c)
{
    assert(data.empty());  // in order that resize fills the vector with zeros.
    data.resize(num_r_c * num_r_c);
    return {data.data(), num_r_c, num_r_c};
}

Eigen::Map<Eigen::VectorXd> toVector(std::vector<double>& data,
                                     std::size_t num_r_c)
{
    assert(data.empty());  // in order that resize fills the vector with zeros.
    data.resize(num_r_c);
    return {data.data(), num_r_c};
}

Eigen::Map<const Eigen::VectorXd> toConstVector(std::vector<double> const& data,
                                                std::size_t num_r_c)
{
    assert(data.size() == num_r_c);
    return {data.data(), num_r_c};
}

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
    auto const num_r_c = local_x_data.size();

    std::vector<double> local_x_perturbed_data(local_x_data);
    auto local_Jac = toMatrix(local_Jac_data, num_r_c);
    auto local_x = toConstVector(local_x_data, num_r_c);
    auto local_xdot = toConstVector(local_xdot_data, num_r_c);

    for (std::size_t i = 0; i < num_r_c; ++i)
    {
        local_x_perturbed_data[i] += eps;
        local_assembler.assemble(t, local_x_perturbed_data, local_M_data,
                                 local_K_data, local_b_data);

        local_x_perturbed_data[i] = local_x[i] - eps;
        local_assembler.assemble(t, local_x_perturbed_data, _local_M_data,
                                 _local_K_data, _local_b_data);

        local_x_perturbed_data[i] = local_x[i];

        auto local_M_p = toMatrix(local_M_data, num_r_c);
        auto local_K_p = toMatrix(local_K_data, num_r_c);
        auto local_b_p = toVector(local_b_data, num_r_c);

        auto local_M_m = toMatrix(_local_M_data, num_r_c);
        auto local_K_m = toMatrix(_local_K_data, num_r_c);
        auto local_b_m = toVector(_local_b_data, num_r_c);

        // TODO check sizes!
        Eigen::MatrixXd const dM_dxi = (local_M_p - local_M_m) / (2.0 * eps);
        Eigen::MatrixXd const dK_dxi = (local_K_p - local_K_m) / (2.0 * eps);
        Eigen::VectorXd const db_dxi = (local_b_p - local_b_m) / (2.0 * eps);
        Eigen::VectorXd const dres_dxi_part =
            dM_dxi * local_xdot + dK_dxi * local_x - db_dxi;

        local_Jac.block(0, i, num_r_c, 1).noalias() = dres_dxi_part;

        local_M_data.clear();
        local_K_data.clear();
        local_b_data.clear();
        _local_M_data.clear();
        _local_K_data.clear();
        _local_b_data.clear();
    }

    local_assembler.assemble(t, local_x_data, local_M_data, local_K_data,
                             local_b_data);

    // TODO no clear!
    auto local_M = toMatrix(local_M_data, num_r_c);
    auto local_K = toMatrix(local_K_data, num_r_c);

    local_Jac += local_M * dxdot_dx + local_K * dx_dx;
}

}  // ProcessLib
