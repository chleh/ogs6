/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_ABSTRACTJACOBIANASSEMBLER_H
#define PROCESSLIB_ABSTRACTJACOBIANASSEMBLER_H

#include "LocalAssemblerInterface.h"

namespace BaseLib
{
class ConfigTree;
}

namespace ProcessLib
{
class AbstractJacobianAssembler
{
public:
    virtual void assembleWithJacobian(std::size_t const mesh_item_id,
                                      LocalAssemblerInterface& local_assembler,
                                      NumLib::LocalToGlobalIndexMap const& dof_table,
                                      const double t, GlobalVector const& x,
                                      GlobalVector const& xdot,
                                      const double dxdot_dx, const double dx_dx,
                                      GlobalMatrix& M, GlobalMatrix& K,
                                      GlobalVector& b,
                                      GlobalMatrix& Jac) const = 0;

    virtual ~AbstractJacobianAssembler() = default;
};

}  // ProcessLib

#endif  // PROCESSLIB_ABSTRACTJACOBIANASSEMBLER_H
