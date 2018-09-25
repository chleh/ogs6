/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "VectorMatrixAssemblerDUNE.h"

#include <cassert>
#include <unordered_map>

#include "LocalAssemblerInterface.h"
#include "NumLib/DOF/DOFTableUtil.h"

#include "Process.h"

namespace ProcessLib
{
VectorMatrixAssemblerDUNE::VectorMatrixAssemblerDUNE(
    std::unique_ptr<AbstractJacobianAssembler>&& jacobian_assembler)
    : _jacobian_assembler(std::move(jacobian_assembler))
{
}

void VectorMatrixAssemblerDUNE::assemble(
    const std::size_t mesh_item_id, LocalAssemblerInterface& local_assembler,
    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>> const&
        dof_tables,
    const double t, const GlobalVector& x, GlobalMatrix& M, GlobalMatrix& K,
    GlobalVector& b, CoupledSolutionsForStaggeredScheme const* const cpl_xs)
{
    std::vector<std::vector<GlobalIndexType>> indices_of_processes;
    indices_of_processes.reserve(dof_tables.size());
    for (std::size_t i = 0; i < dof_tables.size(); i++)
    {
        indices_of_processes.emplace_back(
            NumLib::getIndices(mesh_item_id, dof_tables[i].get()));
    }

    auto const& indices = (cpl_xs == nullptr)
                              ? indices_of_processes[0]
                              : indices_of_processes[cpl_xs->process_id];

    if (cpl_xs == nullptr)
    {
        auto const local_x = x.get(indices);
        local_assembler.assemble(t, local_x, _local_M_data, _local_K_data,
                                 _local_b_data);
    }
    else
    {
        auto local_coupled_xs0 =
            getPreviousLocalSolutions(*cpl_xs, indices_of_processes);
        auto local_coupled_xs =
            getCurrentLocalSolutions(*cpl_xs, indices_of_processes);

        ProcessLib::LocalCoupledSolutions local_coupled_solutions(
            cpl_xs->dt, cpl_xs->process_id, std::move(local_coupled_xs0),
            std::move(local_coupled_xs));

        local_assembler.assembleForStaggeredScheme(t, _local_M_data,
                                                   _local_K_data, _local_b_data,
                                                   local_coupled_solutions);
    }

    auto const num_r_c = indices.size();
    auto const r_c_indices =
        NumLib::LocalToGlobalIndexMap::RowColumnIndices(indices, indices);

    if (!_local_M_data.empty())
    {
        auto const local_M = MathLib::toMatrix(_local_M_data, num_r_c, num_r_c);
        M.add(r_c_indices, local_M);
    }
    if (!_local_K_data.empty())
    {
        auto const local_K = MathLib::toMatrix(_local_K_data, num_r_c, num_r_c);
        K.add(r_c_indices, local_K);
    }
    if (!_local_b_data.empty())
    {
        assert(_local_b_data.size() == num_r_c);
        b.add(indices, _local_b_data);
    }
}

}  // namespace ProcessLib
