/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LocalAssemblerInterface.h"
#include <cassert>
#include "NumLib/DOF/DOFTableUtil.h"

#include "CoupledSolutionsForStaggeredScheme.h"

namespace ProcessLib
{
void LocalAssemblerInterface::assemble(double const /*t*/,
                                       std::vector<double> const& /*local_x*/,
                                       std::vector<double>& /*local_M_data*/,
                                       std::vector<double>& /*local_K_data*/,
                                       std::vector<double>& /*local_b_data*/)
{
    OGS_FATAL(
        "The assemble() function is not implemented in the local assembler.");
}

void LocalAssemblerInterface::assembleForStaggeredScheme(
    double const /*t*/,
    std::vector<double>& /*local_M_data*/,
    std::vector<double>& /*local_K_data*/,
    std::vector<double>& /*local_b_data*/,
    LocalCoupledSolutions const& /*coupled_solutions*/)
{
    OGS_FATAL(
        "The assembleForStaggeredScheme() function is not implemented in the "
        "local assembler.");
}

void LocalAssemblerInterface::assembleWithJacobian(
    double const /*t*/, std::vector<double> const& /*local_x*/,
    std::vector<double> const& /*local_xdot*/, const double /*dxdot_dx*/,
    const double /*dx_dx*/, std::vector<double>& /*local_M_data*/,
    std::vector<double>& /*local_K_data*/,
    std::vector<double>& /*local_b_data*/,
    std::vector<double>& /*local_Jac_data*/)
{
    OGS_FATAL(
        "The assembleWithJacobian() function is not implemented in the local "
        "assembler.");
}

void LocalAssemblerInterface::assembleWithJacobianForStaggeredScheme(
    double const /*t*/, std::vector<double> const& /*local_xdot*/,
    const double /*dxdot_dx*/, const double /*dx_dx*/,
    std::vector<double>& /*local_M_data*/,
    std::vector<double>& /*local_K_data*/,
    std::vector<double>& /*local_b_data*/,
    std::vector<double>& /*local_Jac_data*/,
    LocalCoupledSolutions const& /*local_coupled_solutions*/)
{
    OGS_FATAL(
        "The assembleWithJacobianForStaggeredScheme() function is not implemented in"
        " the local assembler.");
}

void LocalAssemblerInterface::computeSecondaryVariable(
    std::size_t const /*mesh_item_id*/,
    NumLib::AbstractDOFTable const& /*dof_table*/, double const /*t*/,
    GlobalVector const& /*x*/,
    CoupledSolutionsForStaggeredScheme const* /*coupled_xs*/)
{
// TODO [DUNE] re-enable
#if 0
    auto const indices = NumLib::getIndices(mesh_item_id, dof_table);

    if (coupled_xs != nullptr)
        return;

    auto const local_x = x.get(indices);
    computeSecondaryVariableConcrete(t, local_x);
#endif
}

void LocalAssemblerInterface::preTimestep(
    std::size_t const /*mesh_item_id*/,
    NumLib::AbstractDOFTable const& /*dof_table*/, GlobalVector const& /*x*/,
    double const /*t*/, double const /*delta_t*/)
{
// TODO [DUNE] re-enable
#if 0
    auto const indices = NumLib::getIndices(mesh_item_id, dof_table);
    auto const local_x = x.get(indices);

    preTimestepConcrete(local_x, t, delta_t);
#endif
}

void LocalAssemblerInterface::postTimestep(
    std::size_t const /*mesh_item_id*/,
    NumLib::AbstractDOFTable const& /*dof_table*/,
    GlobalVector const& /*x*/)
{
// TODO [DUNE] re-enable
#if 0
    auto const indices = NumLib::getIndices(mesh_item_id, dof_table);
    auto const local_x = x.get(indices);

    postTimestepConcrete(local_x);
#endif
}

void LocalAssemblerInterface::postNonLinearSolver(
    std::size_t const /*mesh_item_id*/,
    NumLib::AbstractDOFTable const& /*dof_table*/, GlobalVector const& /*x*/,
    double const /*t*/, bool const /*use_monolithic_scheme*/)
{
    // TODO [DUNE] re-enable
#if 0
    auto const indices = NumLib::getIndices(mesh_item_id, dof_table);
    auto const local_x = x.get(indices);

    postNonLinearSolverConcrete(local_x, t, use_monolithic_scheme);
#endif
}

}  // namespace ProcessLib
