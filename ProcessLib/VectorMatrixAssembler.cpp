/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "VectorMatrixAssembler.h"

#include <cassert>

#include "LocalAssemblerInterface.h"
#include "NumLib/DOF/DOFTableUtil.h"

namespace ProcessLib
{
void VectorMatrixAssembler::assemble(
    const std::size_t mesh_item_id, LocalAssemblerInterface& local_assembler,
    const NumLib::LocalToGlobalIndexMap& dof_table, const double t,
    const GlobalVector& x, GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    auto const indices = NumLib::getIndices(mesh_item_id, dof_table);
    auto const local_x = x.get(indices);

    _local_M_data.clear();
    _local_K_data.clear();
    _local_b_data.clear();

    local_assembler.assemble(t, local_x, _local_M_data, _local_K_data, _local_b_data);

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

void VectorMatrixAssembler::assembleWithJacobian(
    std::size_t const mesh_item_id, LocalAssemblerInterface& local_assembler,
    NumLib::LocalToGlobalIndexMap const& dof_table, const double t,
    GlobalVector const& x, GlobalVector const& xdot, const double dxdot_dx,
    const double dx_dx, GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b,
    GlobalMatrix& Jac)
{
    auto const indices = NumLib::getIndices(mesh_item_id, dof_table);
    auto const local_x = x.get(indices);

    local_assembler.assembleJacobian(t, local_x, _local_Jac_data);

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

}  // ProcessLib
