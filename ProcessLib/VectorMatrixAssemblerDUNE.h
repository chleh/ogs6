/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>
#include "AbstractJacobianAssembler.h"
#include "BaseLib/Error.h"
#include "CoupledSolutionsForStaggeredScheme.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/NumericsConfig.h"

namespace NumLib
{
class LocalToGlobalIndexMap;
}  // namespace NumLib

namespace ProcessLib
{
class LocalAssemblerInterface;

//! Utility class used to assemble global matrices and vectors.
//!
//! The methods of this class get the global matrices and vectors as input and
//! pass only local data on to the local assemblers.
class VectorMatrixAssemblerDUNE final
{
public:
    explicit VectorMatrixAssemblerDUNE(
        std::unique_ptr<AbstractJacobianAssembler>&& jacobian_assembler);

    //! Assembles\c M, \c K, and \c b.
    //! \remark Jacobian is not assembled here, see assembleWithJacobian().
    void assemble(
        std::size_t const mesh_item_id,
        LocalAssemblerInterface& local_assembler,
        const std::vector<
            std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>& dof_tables,
        double const t, GlobalVector const& x, GlobalMatrix& M, GlobalMatrix& K,
        GlobalVector& b,
        CoupledSolutionsForStaggeredScheme const* const cpl_xs);

    //! Assembles \c M, \c K, \c b, and the Jacobian \c Jac of the residual.
    //! \note The Jacobian must be assembled.
    template <typename Element, typename Basis, typename LocalView,
              typename LocalIndexSet>
    void assembleWithJacobian(
        Element const& e, LocalAssemblerInterface& local_assembler,
        Basis const& basis, LocalView const& localView,
        LocalIndexSet const& localIndexSet, const double t,
        GlobalVector const& x, GlobalVector const& xdot, const double dxdot_dx,
        const double dx_dx, GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b,
        GlobalMatrix& Jac,
        CoupledSolutionsForStaggeredScheme const* const cpl_xs);

private:
    // temporary data only stored here in order to avoid frequent memory
    // reallocations.
    std::vector<double> _local_M_data;
    std::vector<double> _local_K_data;
    std::vector<double> _local_b_data;
    std::vector<double> _local_Jac_data;

    //! Used to assemble the Jacobian.
    std::unique_ptr<AbstractJacobianAssembler> _jacobian_assembler;
};

template <typename Element, typename Basis, typename LocalView,
          typename LocalIndexSet>
void VectorMatrixAssemblerDUNE::assembleWithJacobian(
    Element const& /*e*/, LocalAssemblerInterface& local_assembler,
    Basis const& /*basis*/, LocalView const& localView,
    LocalIndexSet const& localIndexSet, const double t, GlobalVector const& x,
    GlobalVector const& xdot, const double dxdot_dx, const double dx_dx,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac,
    CoupledSolutionsForStaggeredScheme const* const cpl_xs)
{
    std::vector<GlobalIndexType> indices;
    indices.reserve(localView.size());
    for (auto i = decltype(localView.size()){0}; i < localView.size(); ++i)
    {
        assert(localIndexSet.index(i).size() == 1);
        indices.push_back(
            static_cast<GlobalIndexType>(localIndexSet.index(i)[0]));
    }

    auto const local_x = x.get(indices);
    auto const local_xdot = xdot.get(indices);

    _local_M_data.clear();
    _local_K_data.clear();
    _local_b_data.clear();
    _local_Jac_data.clear();

    if (!cpl_xs)
    {
        _jacobian_assembler->assembleWithJacobian(
            local_assembler, t, local_x, local_xdot, dxdot_dx, dx_dx,
            _local_M_data, _local_K_data, _local_b_data, _local_Jac_data);
    }
    else
    {
        // TODO [DUNE] re-enable
#if 0
        auto local_coupled_xs0 = getPreviousLocalSolutionsOfCoupledProcesses(
            *coupling_term, indices);
        auto local_coupled_xs = getCurrentLocalSolutionsOfCoupledProcesses(
            coupling_term->coupled_xs, indices);
        if (local_coupled_xs0.empty() || local_coupled_xs.empty())
        {
            _jacobian_assembler->assembleWithJacobian(
                local_assembler, t, local_x, local_xdot, dxdot_dx, dx_dx,
                _local_M_data, _local_K_data, _local_b_data, _local_Jac_data);
        }
        else
        {
            ProcessLib::LocalCouplingTerm local_coupling_term(
                coupling_term->dt, coupling_term->coupled_processes,
                std::move(local_coupled_xs0), std::move(local_coupled_xs));

            _jacobian_assembler->assembleWithJacobianAndCoupling(
                local_assembler, t, local_x, local_xdot, dxdot_dx, dx_dx,
                _local_M_data, _local_K_data, _local_b_data, _local_Jac_data,
                local_coupling_term);
        }
#endif
        OGS_FATAL("Staggered couppling not implemented, yet.");
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
    if (!_local_Jac_data.empty())
    {
        auto const local_Jac =
            MathLib::toMatrix(_local_Jac_data, num_r_c, num_r_c);
        Jac.add(r_c_indices, local_Jac);
    }
    else
    {
        OGS_FATAL(
            "No Jacobian has been assembled! This might be due to programming "
            "errors in the local assembler of the current process.");
    }
}

}  // namespace ProcessLib
