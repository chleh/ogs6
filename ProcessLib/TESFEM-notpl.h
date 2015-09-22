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

#include <memory>
#include <Eigen/Eigen>

#include "TESProcess-notpl.h"

#include "MathLib/ODE/OdeSolver.h"


namespace ProcessLib
{

namespace TES
{

enum class SecondaryVariables { SOLID_DENSITY, REACTION_RATE,
                                VELOCITY_X, VELOCITY_Y, VELOCITY_Z,
                                REACTION_KINETIC_INDICATOR };

const double NONLINEAR_ERROR_TOLERANCE = 1e-6;



class LADataNoTpl
{
public:
    typedef Eigen::Ref<const Eigen::MatrixXd> MatRef;
    typedef Eigen::Ref<const Eigen::VectorXd> VecRef;

    void assembleIntegrationPoint(
            unsigned integration_point,
            Eigen::MatrixXd* localA,
            Eigen::VectorXd* localRhs,
            std::vector<double> const& localX,
            VecRef const& smN,
            MatRef const& smDNdx,
            const double smDetJ,
            const double weight
            );

    // TESProcessInterface const* _process;
    AssemblyParams const* _AP;

    void init(const unsigned num_int_pts, const unsigned dimension);

    void preEachAssemble();
    void postEachAssemble(Eigen::MatrixXd* localA, Eigen::VectorXd* localRhs,
                          const Eigen::VectorXd& oldX);

    std::vector<double> const&
    getIntegrationPointValues(SecondaryVariables var) const;

private:
    Eigen::Matrix3d getMassCoeffMatrix(const unsigned int_pt);
    Eigen::MatrixXd getLaplaceCoeffMatrix(const unsigned int_pt, const unsigned dim);
    Eigen::Matrix3d getAdvectionCoeffMatrix(const unsigned int_pt);
    Eigen::Matrix3d getContentCoeffMatrix(const unsigned int_pt);
    Eigen::Vector3d getRHSCoeffVector(const unsigned int_pt);

    void preEachAssembleIntegrationPoint(
            const unsigned int_pt,
            std::vector<double> const& localX,
            VecRef const& smN,
            MatRef const& smDNdx
            );

    void initReaction(
            const unsigned int_pt,
            std::vector<double> const& localX);

    // many values taken from zeolite-adsorption-benchmark-snap/start-at-0.99


    // nodal quantities, secondary variables
    std::vector<double> _solid_density;
    std::vector<double> _solid_density_prev_ts;

    std::vector<double> _reaction_rate; // dC/dt * _rho_SR_dry
    std::vector<double> _reaction_rate_prev_ts; // could also be calculated from _solid_density_prev_ts

    std::vector<std::vector<double> > _velocity;
    // std::vector<double> _velocity_x;
    // std::vector<double> _velocity_x;
    // Eigen::MatrixXd _velocity; // row index: gauss point, column index: dimension x/y/z

    std::vector<double> _reaction_rate_indicator;

    // integration point values of unknowns
    double _p = -888.888; // gas pressure
    double _T = -888.888; // temperature
    double _vapour_mass_fraction = -888.888;     // fluid mass fraction of the second component

    // temporary storage for some properties
    // values do not change during the assembly of one integration point
    double _rho_GR = -888.888;
    double _p_V = -888.888; // vapour partial pressure
    double _qR = 88888.88888;  // reaction rate, use this in assembly!!!

    std::unique_ptr<Eigen::MatrixXd> _Lap;
    std::unique_ptr<Eigen::MatrixXd> _Mas;
    std::unique_ptr<Eigen::MatrixXd> _Adv;
    std::unique_ptr<Eigen::MatrixXd> _Cnt;
    std::unique_ptr<Eigen::VectorXd> _rhs;

    std::unique_ptr<MathLib::OdeSolver<2> > _ode_solver;
};


void
ogs5OutVec(const LADataNoTpl::VecRef& vec);

void
ogs5OutMat(const LADataNoTpl::MatRef& vec);


} // namespace TES

} // namespace ProcessLib

#endif // PROCESS_LIB_TES_FEM_NOTPL_H_
