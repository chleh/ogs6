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

#include <iostream>

#include "TESFEM-notpl.h"


const double GAS_CONST = 8.3144621;

template <typename T>
T square(const T& v)
{
	return v * v;
}



static inline double fluid_cp(const double p, const double T, const double x)
{
	return 888.888 + 0.0 * (p + T + x);
}

static inline double fluid_rho(const double p, const double T, const double x)
{
	return 888.888 + 0.0 * (p + T + x);
}

static inline double fluid_eta(const double p, const double T, const double x)
{
	return 888.888 + 0.0 * (p + T + x);
}

static inline double fluid_heat_cond(const double p, const double T, const double x)
{
	return 888.888 + 0.0 * (p + T + x);
}

static inline double solid_cp(const double rho_SR)
{
	return 888.888 + 0.0 * (rho_SR);
}

static inline double evap_enthalp(const double p, const double T, const double x)
{
	return 888.888 + 0.0 * (p + T + x);
}


namespace ProcessLib
{

namespace TES
{

Eigen::Matrix3d
LADataNoTpl::
getMassCoeffMatrix()
{
	const double rho_GR = fluid_rho(_p, _T, _x);

	const double cpG   = fluid_cp(_p, _T, _x);
	const double cpS   = solid_cp(_rho_SR);

	double dxn_dxm = _M_inert * _M_react // 0 is inert, 1 is reactive
					 / square(_M_inert * _x + _M_react * (1.0 - _x));

	const double M_pp = _poro/_p * rho_GR;
	const double M_pT = -_poro/_T *  rho_GR;
	const double M_px = (_M_react-_M_inert) * _p / (GAS_CONST * _T) * dxn_dxm * _poro;

	const double M_Tp = -_poro;
	const double M_TT = _poro * rho_GR * cpG + (1.0-_poro) * _rho_SR * cpS;
	const double M_Tx = 0.0;

	const double M_xp = 0.0;
	const double M_xT = 0.0;
	const double M_xx = _poro * rho_GR;


	Eigen::Matrix3d M;
	M << M_pp, M_pT, M_px,
		 M_Tp, M_TT, M_Tx,
		 M_xp, M_xT, M_xx;

	return M;
}


Eigen::MatrixXd
LADataNoTpl::
getLaplaceCoeffMatrix(const unsigned dim)
{
	// TODO implement

	const double rho_GR = fluid_rho(_p, _T, _x);
	const double eta_GR = fluid_eta(_p, _T, _x);

	const double lambda_F = fluid_heat_cond(_p, _T, _x);

	// TODO: k_rel
	auto const L_pp = _solid_perm_tensor.block(0,0,dim,dim) * rho_GR / eta_GR;

	auto const L_pT = Eigen::MatrixXd::Zero(dim, dim);
	auto const L_px = Eigen::MatrixXd::Zero(dim, dim);

	auto const L_Tp = Eigen::MatrixXd::Zero(dim, dim);

	// TODO: add zeolite part
	auto const L_TT = Eigen::MatrixXd::Identity(dim, dim)
					  * ( _poro * lambda_F + (1.0 - _poro) * _solid_heat_cond);

	auto const L_Tx = Eigen::MatrixXd::Zero(dim, dim);

	auto const L_xp = Eigen::MatrixXd::Zero(dim, dim);
	auto const L_xT = Eigen::MatrixXd::Zero(dim, dim);

	auto const L_xx = Eigen::MatrixXd::Identity(dim, dim)
					  * _tortuosity * _poro * rho_GR * _diffusion_coefficient_component;


	Eigen::MatrixXd L(dim*3, dim*3);

	L.block(0, 0, dim, dim) = L_pp;
	L.block(0, 1, dim, dim) = L_pT;
	L.block(0, 2, dim, dim) = L_px;

	L.block(1, 0, dim, dim) = L_Tp;
	L.block(1, 1, dim, dim) = L_TT;
	L.block(1, 2, dim, dim) = L_Tx;

	L.block(2, 0, dim, dim) = L_xp;
	L.block(2, 1, dim, dim) = L_xT;
	L.block(2, 2, dim, dim) = L_xx;

	return L;
}


Eigen::Matrix3d
LADataNoTpl::
getAdvectionCoeffMatrix()
{
	const double rho_GR = fluid_rho(_p, _T, _x);

	const double cpG   = fluid_cp(_p, _T, _x);

	const double A_pp = 0.0;
	const double A_pT = 0.0;

	const double A_px = 0.0;

	const double A_Tp = 0.0;

	const double A_TT = rho_GR * cpG; // porosity?
	const double A_Tx = 0.0;

	const double A_xp = 0.0;
	const double A_xT = 0.0;
	const double A_xx = rho_GR; // porosity?


	Eigen::Matrix3d A;
	A << A_pp, A_pT, A_px,
		 A_Tp, A_TT, A_Tx,
		 A_xp, A_xT, A_xx;

	return A;
}


Eigen::Matrix3d
LADataNoTpl::
getContentCoeffMatrix()
{
	const double C_pp = 0.0;
	const double C_pT = 0.0;

	const double C_px = 0.0;

	const double C_Tp = 0.0;

	const double C_TT = 0.0;
	const double C_Tx = 0.0;

	const double C_xp = 0.0;
	const double C_xT = 0.0;
	const double C_xx = (1.0 - _poro) * _reaction_rate;


	Eigen::Matrix3d C;
	C << C_pp, C_pT, C_px,
		 C_Tp, C_TT, C_Tx,
		 C_xp, C_xT, C_xx;

	return C;
}


Eigen::Vector3d
LADataNoTpl::
getRHSCoeffVector()
{
	const double rho_GR = fluid_rho(_p, _T, _x);

	const double rhs_p = (_poro - 1.0) * _reaction_rate;

	const double rhs_T = rho_GR * _fluid_specific_heat_source
						 + (1.0 - _poro) * _reaction_rate * evap_enthalp(_p, _T, _x)
						 + _rho_SR * (1.0 - _poro) * _solid_specific_heat_source;

	const double rhs_x = (_poro - 1.0) * _reaction_rate; // TODO: what if x < 0.0


	Eigen::Vector3d rhs;
	rhs << rhs_p,
		 rhs_T,
		 rhs_x;

	return rhs;
}


void
LADataNoTpl::
init(unsigned GlobalDim)
{
	_GlobalDim = GlobalDim;
}


void
LADataNoTpl::
assembleIntegrationPoint(
		Eigen::MatrixXd* localA, Eigen::VectorXd* localRhs,
		const MatRef &smN, const MatRef &smDNdx, const double smDetJ,
		const double weight)
{
    std::cerr << "localA:\n" << (*localA) << std::endl;
    std::cerr << "localA block\n" << (localA->template block<2,2>(0,0)) << std::endl;
    std::cerr << "coeff:\n" << smDNdx.transpose() * 1.0 * smDNdx * smDetJ * weight << std::endl;

    auto const N = smN.size();
    auto const D = _GlobalDim;
    unsigned const NODAL_DOF = 3;

    auto const laplaceCoeffMat = getLaplaceCoeffMatrix(D);
    auto const massCoeffMat    = getMassCoeffMatrix();
    auto const advCoeffMat     = getAdvectionCoeffMatrix();
    auto const contentCoeffMat = getContentCoeffMatrix();

    Eigen::MatrixXd const velocity = Eigen::MatrixXd::Constant(D, 1, 888.888);

    for (unsigned r=0; r<NODAL_DOF; ++r)
    {
        for (unsigned c=0; c<NODAL_DOF; ++c)
        {
            std::cerr << "vel prod: " << smDNdx.transpose() * velocity << std::endl;
            localA->block(N*r, N*c, N, N).noalias() +=
                    (
                        smDNdx.transpose() * laplaceCoeffMat.block(D*r, D*c, D, D) * smDNdx
                        + smN * (massCoeffMat(r, c) + contentCoeffMat(r, c)) * smN.transpose()
                        + smN * advCoeffMat(r, c) * velocity.transpose() * smDNdx
                    )
                    * smDetJ * weight;

            std::cerr << "localA block(" << N*r << "," << N*c << ")\n"
                      << (localA->block(N*r, N*c, N, N)) << std::endl;
            std::cerr << "coeff:\n"
                      << smN.transpose() * massCoeffMat(r, c) * smN
                         * smDetJ * weight << std::endl;
        }
    }

    std::cerr << "sm.N: %s" << smN << std::endl;
    std::cerr << "sm.dNdx: %s" << smDNdx << std::endl;
    std::cerr << "mass coeffs:\n" << massCoeffMat << std::endl;

    auto const rhsCoeffVector = getRHSCoeffVector();

    for (unsigned r=0; r<NODAL_DOF; ++r)
    {
        localRhs->block(N*r, 0, N, 1) +=
                rhsCoeffVector(r) * smN * smDetJ * weight;
    }
}

} // namespace TES

} // namespace ProcessLib

