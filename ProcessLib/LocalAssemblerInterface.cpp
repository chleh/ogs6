/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <cassert>
#include "LocalAssemblerInterface.h"
#include "NumLib/DOF/DOFTableUtil.h"

namespace ProcessLib
{
void LocalAssemblerInterface::assemble(
    const std::size_t mesh_item_id,
    const NumLib::LocalToGlobalIndexMap& dof_table, const double t,
    const GlobalVector& x, GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    auto const indices = NumLib::getIndices(mesh_item_id, dof_table);
    auto const local_x = x.get(indices);

    _local_M_data.clear();
    _local_K_data.clear();
    _local_b_data.clear();

    assembleConcrete(t, local_x, _local_M_data, _local_K_data, _local_b_data);

    auto const num_r_c = indices.size();
    auto const r_c_indices =
        NumLib::LocalToGlobalIndexMap::RowColumnIndices(indices, indices);

    if (!_local_M_data.empty()) {
        assert(_local_M_data.size() == num_r_c * num_r_c);
        auto const local_M =
            Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                                     Eigen::RowMajor>>(_local_M_data.data(),
                                                       num_r_c, num_r_c);
        M.add(r_c_indices, local_M);
    }
    if (!_local_K_data.empty()) {
        assert(_local_K_data.size() == num_r_c * num_r_c);
        auto const local_K =
            Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                                     Eigen::RowMajor>>(_local_K_data.data(),
                                                       num_r_c, num_r_c);
        K.add(r_c_indices, local_K);
    }
    if (!_local_b_data.empty()) {
        assert(_local_b_data.size() == num_r_c);
        auto const local_b =
            Eigen::Map<Eigen::VectorXd>(_local_b_data.data(), num_r_c);
        b.add(indices, local_b);
    }
}

void LocalAssemblerInterface::assembleJacobian(
    const std::size_t mesh_item_id,
    const NumLib::LocalToGlobalIndexMap& dof_table, const double t,
    const GlobalVector& x, GlobalMatrix& Jac)
{
    auto const indices = NumLib::getIndices(mesh_item_id, dof_table);
    auto const local_x = x.get(indices);

    assembleJacobianConcrete(t, local_x, _local_Jac_data);

    if (!_local_Jac_data.empty()) {
        auto const num_r_c = indices.size();
        auto const r_c_indices =
            NumLib::LocalToGlobalIndexMap::RowColumnIndices(indices, indices);
        assert(_local_Jac_data.size() == num_r_c * num_r_c);
        auto const local_Jac = Eigen::Map<Eigen::MatrixXd>(
            _local_Jac_data.data(), num_r_c, num_r_c);
        Jac.add(r_c_indices, local_Jac);
    }
}

void LocalAssemblerInterface::assembleJacobianConcrete(
        double const /*t*/, std::vector<double> const& /*local_x*/,
        std::vector<double>& /*local_Jac_data*/)
{
    OGS_FATAL(
        "The assembleJacobian() function is not implemented in the local "
        "assembler.");
}

void LocalAssemblerInterface::preTimestep(
    std::size_t const mesh_item_id,
    NumLib::LocalToGlobalIndexMap const& dof_table, GlobalVector const& x,
    double const t, double const delta_t)
{
    auto const indices = NumLib::getIndices(mesh_item_id, dof_table);
    auto const local_x = x.get(indices);

    preTimestepConcrete(local_x, t, delta_t);
}

void LocalAssemblerInterface::postTimestep(
    std::size_t const mesh_item_id,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    GlobalVector const& x)
{
    auto const indices = NumLib::getIndices(mesh_item_id, dof_table);
    auto const local_x = x.get(indices);

    postTimestepConcrete(local_x);
}

}  // namespace ProcessLib
