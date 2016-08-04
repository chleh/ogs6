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
class AnalyticalJacobianAssembler : public AbstractJacobianAssembler
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
        // local_assembler.assemble(mesh_item_id, dof_table, t, x, M, K, b);
        // local_assembler.assembleJacobian(mesh_item_id, dof_table, t, x, Jac);
    }
};

}  // ProcessLib

#endif  // PROCESSLIB_ANALYTICALJACOBIANASSEMBLER_H
