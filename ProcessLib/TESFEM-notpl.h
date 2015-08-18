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

#include "TESProcess-notpl.h"


namespace ProcessLib
{

namespace TES
{

const unsigned NODAL_DOF = 3;
const unsigned NODAL_DOF_2ND = 2; // loading or solid density, and reaction rate


const double M_N2  = 0.028013;
const double M_H2O = 0.018016;


class LADataNoTpl
{
public:
    typedef Eigen::Ref<const Eigen::MatrixXd> MatRef;
    typedef Eigen::Ref<const Eigen::VectorXd> VecRef;

    void assembleIntegrationPoint(
            Eigen::MatrixXd* localA,
            Eigen::VectorXd* localRhs,
            std::vector<double> const& localX,
            std::vector<double> const& localSecondaryVariables,
            VecRef const& smN,
            MatRef const& smDNdx,
            const double smDetJ,
            const double weight
            );


    TESProcessInterface const* _process;


private:
    Eigen::Matrix3d getMassCoeffMatrix();
    Eigen::MatrixXd getLaplaceCoeffMatrix(const unsigned dim);
    Eigen::Matrix3d getAdvectionCoeffMatrix();
    Eigen::Matrix3d getContentCoeffMatrix();
    Eigen::Vector3d getRHSCoeffVector();

    void preEachAssembleIntegrationPoint(
            std::vector<double> const& localX,
            std::vector<double> const& localSecondaryVariables,
            VecRef const& smN,
            MatRef const& smDNdx
            );

    // many values taken from zeolite-adsorption-benchmark-snap/start-at-0.99

    double _fluid_specific_heat_source = 0.0;
    double _cpG = 1012.0; // specific isobaric fluid heat capacity

    Eigen::MatrixXd _solid_perm_tensor = Eigen::MatrixXd::Identity(3, 3) * 6.94e-14; // TODO get dimensions
    double _solid_specific_heat_source = 0.0;
    double _solid_heat_cond = 0.4;
    double _cpS = 880.0;    // specific isobaric solid heat capacity

    double _tortuosity = 1.0;
    double _diffusion_coefficient_component = 9.65e-5; // ???

    double _poro = 0.7;

    const double _rho_SR_dry = 1150.0;

    double _M_inert = M_N2; // N2
    double _M_react = M_H2O;

    // integration point values of unknowns
    double _p = 888.888; // gas pressure
    double _T = 888.888; // temperature
    double _vapour_mass_fraction = 0.5;     // fluid mass fraction of the second component

    // integration point values of secondary veriables
    double _solid_density = 888.888; // rho_SR
    double _reaction_rate = 888.888; // dC/dt * _rho_SR_dry

    // temporary storage for some properties
    // values do not change during the assembly of one integration point
    double _rho_GR = 888.888;
    double _p_V = 888.888; // vapour partial pressure
    // double _vapour_molar_fraction = 888.888;
};


} // namespace TES

} // namespace ProcessLib

#endif // PROCESS_LIB_TES_FEM_NOTPL_H_
