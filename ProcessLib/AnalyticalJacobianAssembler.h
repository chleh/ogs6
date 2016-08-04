/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_ANALYTICALJACOBIANASSEMBLER_H
#define PROCESSLIB_ANALYTICALJACOBIANASSEMBLER_H

#include "AbstractJacobianAssembler.h"

namespace BaseLib
{
class ConfigTree;
}

namespace ProcessLib
{
template <typename LocalAssemblerInterface>
class AnalyticalJacobianAssembler
    : public AbstractJacobianAssembler<LocalAssemblerInterface>
{
public:
    void assembleWithJacobian(std::size_t const mesh_item_id,
                              LocalAssemblerInterface& local_assembler,
                              NumLib::LocalToGlobalIndexMap const& dof_table,
                              const double t, GlobalVector const& x,
                              GlobalVector const& xdot, const double dxdot_dx,
                              const double dx_dx, GlobalMatrix& M,
                              GlobalMatrix& K, GlobalVector& b,
                              GlobalMatrix& Jac) const override
    {
        // assemble_jacobian_callback(t, x, xdot, dxdot_dx, dx_dx, M, K, b, Jac);
    }
};

}  // ProcessLib

#endif  // PROCESSLIB_ANALYTICALJACOBIANASSEMBLER_H
