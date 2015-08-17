/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * The code of this file is used to decouple the evaluation of matrix elements from the rest of OGS6,
 * not all of OGS6 has to be recompiled every time a small change is done.
 */

#ifndef PROCESS_LIB_TES_FEM_NOTPL_H_
#define PROCESS_LIB_TES_FEM_NOTPL_H_

#include <Eigen/Eigen>

namespace ProcessLib
{

namespace TES
{

const unsigned NODAL_DOF = 3;


struct LADataNoTpl
{
    typedef Eigen::Ref<const Eigen::MatrixXd> MatRef;
    typedef Eigen::Ref<const Eigen::VectorXd> VecRef;

    Eigen::Matrix3d getMassCoeffMatrix();
    Eigen::MatrixXd getLaplaceCoeffMatrix(const unsigned dim);
    Eigen::Matrix3d getAdvectionCoeffMatrix();
    Eigen::Matrix3d getContentCoeffMatrix();
    Eigen::Vector3d getRHSCoeffVector();

    void assembleIntegrationPoint(
            Eigen::MatrixXd* localA,
            Eigen::VectorXd* localRhs,
            std::vector<double> const& localX,
            VecRef const& smN,
            MatRef const& smDNdx,
            const double smDetJ,
            const double weight
            );

    double _reaction_rate = 888.888;

    double _fluid_specific_heat_source = 888.888;

    double _solid_specific_heat_source = 888.888;
    Eigen::MatrixXd _solid_perm_tensor = Eigen::MatrixXd::Constant(3, 3, 888.888); // TODO get dimensions
    double _solid_heat_cond = 888.888;

    double _tortuosity = 888.888;
    double _diffusion_coefficient_component = 888.888;

    double _poro = 888.888;

    double _p = 888.888;
    double _T = 888.888;
    double _x = 0.5;

    double _rho_SR = 888.888;

    double _M_inert = 888.888;
    double _M_react = 888.888;
};


} // namespace TES

} // namespace ProcessLib

#endif // PROCESS_LIB_TES_FEM_NOTPL_H_
