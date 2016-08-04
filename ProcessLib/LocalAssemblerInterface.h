/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_LOCALASSEMBLERINTERFACE_H
#define PROCESSLIB_LOCALASSEMBLERINTERFACE_H

#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/NumericsConfig.h"

namespace ProcessLib
{

/*! Common interface for local assemblers
 * NumLib::ODESystemTag::FirstOrderImplicitQuasilinear ODE systems.
 *
 * \todo Generalize to other NumLib::ODESystemTag's.
 */
class LocalAssemblerInterface
{
public:
    virtual ~LocalAssemblerInterface() = default;

    void assemble(std::size_t const mesh_item_id,
                  NumLib::LocalToGlobalIndexMap const& dof_table,
                  double const t, GlobalVector const& x, GlobalMatrix& M,
                  GlobalMatrix& K, GlobalVector& b);

    void assembleJacobian(std::size_t const mesh_item_id,
                          NumLib::LocalToGlobalIndexMap const& dof_table,
                          double const t, GlobalVector const& x,
                          GlobalMatrix& Jac);

    virtual void preTimestep(std::size_t const mesh_item_id,
                             NumLib::LocalToGlobalIndexMap const& dof_table,
                             GlobalVector const& x, double const t,
                             double const delta_t);

    virtual void postTimestep(std::size_t const mesh_item_id,
                              NumLib::LocalToGlobalIndexMap const& dof_table,
                              GlobalVector const& x);

protected:
    virtual void assembleConcrete(
        double const t, std::vector<double> const& local_x,
        std::vector<double>& local_M_data, std::vector<double>& local_K_data,
        std::vector<double>& local_b_data) = 0;

    virtual void assembleJacobianConcrete(
        double const t, std::vector<double> const& local_x,
        std::vector<double>& local_Jac_data);

    virtual void preTimestepConcrete(std::vector<double> const& /*local_x*/,
                                     double const /*t*/, double const /*dt*/)
    {
    }

    virtual void postTimestepConcrete(std::vector<double> const& /*local_x*/) {}

private:
    std::vector<double> _local_M_data;
    std::vector<double> _local_K_data;
    std::vector<double> _local_b_data;
    std::vector<double> _local_Jac_data;
};

} // namespace ProcessLib

#endif // PROCESSLIB_LOCALASSEMBLERINTERFACE_H
