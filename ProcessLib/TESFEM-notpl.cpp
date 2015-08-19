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
#include <cstdio>

#include "logog/include/logog.hpp"

#include "TESFEM-notpl.h"


const double GAS_CONST = 8.3144621;




template <typename T>
T square(const T& v)
{
	return v * v;
}


#if 0
static double fluid_specific_isobaric_heat_capacity(
		const double /*p*/, const double /*T*/, const double /*x*/)
{
	// OGS-5 heat cap model 1 value 1000
	// constant cp

	return 1000.0;
}
#endif


static double fluid_density(const double p, const double T, const double x)
{
	// OGS-5 density model 26

	const double M0 = ProcessLib::TES::M_N2;
	const double M1 = ProcessLib::TES::M_H2O;

	const double xn = M0*x/(M0*x + M1*(1.0-x));

	return p / (GAS_CONST * T) * (M1*xn + M0*(1.0-xn));;
}

static double fluid_viscosity_N2(double rho, const double T)
{
	const double rho_c = 314;             // [kg/m3]
	const double CVF = 14.058;            // [1e-3 Pa-s]

	const double sigma = 0.36502496e-09;
	const double k = 1.38062e-23;
	const double eps = 138.08483e-23;
	const double c1 = 0.3125;
	const double c2 = 2.0442e-49;

	double A[5],C[5];

	A[0] = 0.46649;
	A[1] = -0.57015;
	A[2] = 0.19164;
	A[3] = -0.03708;
	A[4] = 0.00241;

	C[0] = -20.09997;
	C[1] = 3.4376416;
	C[2] = -1.4470051;
	C[3] = -0.027766561;
	C[4] = -0.21662362;

	const double T_star = T * k / eps;
	rho = rho / rho_c;

	double Omega = 0.0;
	for (unsigned i = 0; i < 5; i++)
		Omega += A[i] * std::pow(log(T_star),i);

	Omega = exp (Omega);

	//eta in [Pa*s]
	const double eta_0 = c1 * sqrt(c2 * T) / (sigma * sigma * Omega);

	double sum = 0.0;
	for (unsigned i = 2; i < 5; i++)
		sum = sum + C[i] * std::pow(rho, i-1);

	//
	const double eta_r = CVF * 1e-6 * (C[0] / (rho - C[1]) + C[0] / C[1] + sum);

	return eta_0 + eta_r;               // [Pa*s]
}

static double fluid_viscosity_H2O(double rho, double T)
{
	double my,my_0,my_1;
	double H[4],h[6][7];
	double sum1 = 0,sum2 = 0,sum3 = 0;
	int i,j;

	T = T / 647.096;
	rho = rho / 322.0;

	H[0] = 1.67752;
	H[1] = 2.20462;
	H[2] = 0.6366564;
	H[3] = -0.241605;

	for (i = 0; i < 6; i++)
		for (j = 0; j < 7; j++)
			h[i][j] = 0;
	h[0][0] = 0.520094000;
	h[1][0] = 0.085089500;
	h[2][0] = -1.083740000;
	h[3][0] = -0.289555000;
	h[0][1] = 0.222531000;
	h[1][1] = 0.999115000;
	h[2][1] = 1.887970000;
	h[3][1] = 1.266130000;
	h[5][1] = 0.120573000;
	h[0][2] = -0.281378000;
	h[1][2] = -0.906851000;
	h[2][2] = -0.772479000;
	h[3][2] = -0.489837000;
	h[4][2] = -0.257040000;
	h[0][3] = 0.161913000;
	h[1][3] = 0.257399000;
	h[0][4] = -0.032537200;
	h[3][4] = 0.069845200;
	h[4][5] = 0.008721020;
	h[3][6] = -0.004356730;
	h[5][6] = -0.000593264;

	for(i = 0; i < 4; i++)
		sum1 = sum1 + (H[i] / std::pow(T,i));

	my_0 = 100 * sqrt(T) / sum1;

	for(i = 0; i < 6; i++)
	{
		for (j = 0; j < 7; j++)
			sum3 = sum3 + h[i][j] * std::pow(rho - 1,j);
		sum2 = sum2 + std::pow(1 / T - 1,i) * sum3;
		sum3 = 0;
	}

	my_1 = exp(rho * sum2);

	my = (my_0 * my_1) / 1e6;
	return my;
}

static double fluid_viscosity(const double p, const double T, const double x)
{
	// OGS 5 viscosity model 26

	const double M0 = ProcessLib::TES::M_N2;
	const double M1 = ProcessLib::TES::M_H2O;

	//reactive component
	const double x0 = M0*x/(M0*x + M1*(1.0-x)); //mass in mole fraction
	const double V0 = fluid_viscosity_H2O(M1*p/(GAS_CONST*T), T);
	//inert component
	const double x1 = 1.0 - x0;
	const double V1 = fluid_viscosity_N2(M0*p/(GAS_CONST*T), T);

	const double M0_over_M1 (M1/M0); //reactive over inert
	const double V0_over_V1 (V0/V1);

	const double phi_12 =   (1.0 + sqrt(V0_over_V1) * pow(1.0/M0_over_M1, 0.25))
						  * (1.0 + sqrt(V0_over_V1) * pow(1.0/M0_over_M1, 0.25))
						  / pow(8.0*(1.0+M0_over_M1),0.5);
	const double phi_21 = phi_12 * M0_over_M1 / V0_over_V1;

	return V0*x0 / (x0 + x1 * phi_12)
			+ V1*x1 / (x1 + x0 * phi_21);
}

static double fluid_heat_conductivity_N2(double rho, double T)
{
	const double X1 = 0.95185202;
	const double X2 = 1.0205422;

	const double rho_c = 314;             // [kg/m3]
	const double M = 28.013;
	const double k = 1.38062e-23;
	const double eps = 138.08483e-23;
	const double N_A = 6.02213E26;
	const double R = 8.31434;
	const double CCF = 4.173;             //mW/m/K

	const double c1 = 0.3125;
	const double c2 = 2.0442e-49;
	const double sigma = 0.36502496e-09;

	double F;
	double A[5],f[9],C[4];
	double sum = 0,eta_0,c_v0,T_star,Omega = 0;
	double lamda_tr,lamda_in,lamda_r,lamda_0,lamda;

	int i;

	T_star = T * k / eps;
	rho = rho / rho_c;

	A[0] = 0.46649;
	A[1] = -0.57015;
	A[2] = 0.19164;
	A[3] = -0.03708;
	A[4] = 0.00241;

	f[0] = -0.837079888737e3;
	f[1] = 0.37914711487e2;
	f[2] = -0.601737844275;
	f[3] = 0.350418363823e1;
	f[4] = -0.874955653028e-5;
	f[5] = 0.148968607239e-7;
	f[6] = -0.256370354277e-11;
	f[7] = 0.100773735767e1;
	f[8] = 0.335340610e4;

	C[0] = 3.3373542;
	C[1] = 0.37098251;
	C[2] = 0.89913456;
	C[3] = 0.16972505;

	// dilute heat conductivity
	for (i = 0; i < 7; i++)
		sum = sum + f[i] * pow(T,(i - 3));
	const double temp (exp ((f[8] / T)) - 1);
	c_v0 = R * (sum + ((f[7] * (f[8] / T) * (f[8] / T) * (exp((f[8] / T)))) / (temp * temp) - 1));
	sum = 0;

	double cvint;
	cvint = c_v0 * 1000 / N_A;

	// dilute gas viscosity
	for (i = 0; i < 5; i++)
		Omega = Omega + A[i] * std::pow(log(T_star),i);
	Omega = exp (Omega);

	//eta in [Pa*s]
	eta_0 = 1e6 * (c1 * sqrt(c2 * T) / (sigma * sigma * Omega));

	F = eta_0 * k * N_A / (M * 1000);

	lamda_tr = 2.5 * (1.5 - X1);
	lamda_in = X2 * (cvint / k + X1);

	lamda_0 = F * (lamda_tr + lamda_in);
	sum = 0;
	for (i = 0; i < 4; i++)
		sum = sum + C[i] * std::pow(rho,(i + 1));

	lamda_r = sum * CCF;

	lamda = (lamda_0 + lamda_r) / 1000;   //lamda in [W/m/K]

	return lamda;
}

static double fluid_heat_conductivity_H2O(double rho, double T)
{
	double lamda,lamda_0,lamda_1,lamda_2;
	double sum1 = 0;
	double S,Q,dT;
	double a[4],b[3],B[2],d[4],C[6];
	int i;

	T = T / 647.096;
	rho = rho / 317.11;

	a[0] =  0.0102811;
	a[1] =  0.0299621;
	a[2] =  0.0156146;
	a[3] = -0.00422464;

	b[0] = -0.397070;
	b[1] =  0.400302;
	b[2] =  1.060000;

	B[0] = -0.171587;
	B[1] =  2.392190;

	d[0] = 0.0701309;
	d[1] = 0.0118520;
	d[2] = 0.00169937;
	d[3] = -1.0200;

	C[0] = 0.642857;
	C[1] = -4.11717;
	C[2] = -6.17937;
	C[3] = 0.00308976;
	C[4] = 0.0822994;
	C[5] = 10.0932;

	for (i = 0; i < 4; i++)
		sum1 = sum1 + a[i] * std::pow(T,i);

	lamda_0 = sqrt(T) * sum1;
	lamda_1 = b[0] + b[1] * rho + b[2] * exp(B[0] * (rho + B[1]) * (rho + B[1]));

	dT = fabs(T - 1) + C[3];
	Q = 2 + (C[4] / pow(dT,3. / 5.));

	if (T >= 1)
		S = 1 / dT;
	else
		S = C[5] / pow(dT,3. / 5.);

	lamda_2 =
	        (d[0] /
	         std::pow(T, 10) + d[1]) * std::pow(rho,9. / 5.) * exp(C[0] * (1 - pow(rho,14. / 5.)))
	        + d[2]* S* std::pow(rho,Q) * exp((Q / (1. + Q)) * (1 - std::pow(rho,(1. + Q))))
	        + d[3] * exp(C[1] * pow(T,3. / 2.) + C[2] / std::pow(rho,5));

	lamda = (lamda_0 + lamda_1 + lamda_2); // lamda in [W/m/K]

	return lamda;
}

static double fluid_heat_conductivity(const double p, const double T, const double x)
{
	// OGS 5 fluid heat conductivity model 11

	const double M0 = ProcessLib::TES::M_N2;
	const double M1 = ProcessLib::TES::M_H2O;

	// TODO [CL] max() is redundant if the fraction is guaranteed to be between 0 and 1.
	//reactive component
	const double x0 = std::max(M0*x/(M0*x + M1*(1.0-x)), 0.); // convert mass to mole fraction
	const double k0 = fluid_heat_conductivity_H2O(M1*p/(GAS_CONST*T), T);
	//inert component
	const double x1 = 1.0 - x0;
	const double k1 = fluid_heat_conductivity_N2(M0*p/(GAS_CONST*T), T);

	const double M1_over_M2 = M1/M0; //reactive over inert
	const double V1_over_V2 = fluid_viscosity_H2O(M1*p/(GAS_CONST*T), T)
							/ fluid_viscosity_N2 (M0*p/(GAS_CONST*T), T);
	const double L1_over_L2 = V1_over_V2 / M1_over_M2;

	const double phi_12 =   (1.0 + pow(L1_over_L2, 0.5) * pow(M1_over_M2, -0.25))
						  * (1.0 + pow(V1_over_V2, 0.5) * pow(M1_over_M2, -0.25))
						  / pow(8.0 * (1.0 + M1_over_M2), 0.5);
	const double phi_21 = phi_12 * M1_over_M2 / V1_over_V2;

	return k0*x0 / (x0+x1*phi_12) + k1*x1 / (x1+x0*phi_21);
}

#if 0
static double solid_specific_isobaric_heat_capacity(const double /*rho_SR*/)
{
	// OGS 5 heat capacity model 1 (constant) value 620.573
	return 620.573;
}
#endif


namespace ProcessLib
{

namespace TES
{

Eigen::Matrix3d
LADataNoTpl::
getMassCoeffMatrix()
{
	double dxn_dxm = _M_inert * _M_react
					 / square(_M_inert * _vapour_mass_fraction + _M_react * (1.0 - _vapour_mass_fraction));

	const double M_pp = _poro/_p * _rho_GR;
	const double M_pT = -_poro/_T *  _rho_GR;
	const double M_px = (_M_react-_M_inert) * _p / (GAS_CONST * _T) * dxn_dxm * _poro;

	const double M_Tp = -_poro;
	const double M_TT = _poro * _rho_GR * _cpG + (1.0-_poro) * _solid_density * _cpS;
	const double M_Tx = 0.0;

	const double M_xp = 0.0;
	const double M_xT = 0.0;
	const double M_xx = _poro * _rho_GR;


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
	const double eta_GR = fluid_viscosity(_p, _T, _vapour_mass_fraction);

	const double lambda_F = fluid_heat_conductivity(_p, _T, _vapour_mass_fraction);
	const double lambda_S = _solid_heat_cond;

	// TODO: k_rel
	auto const L_pp = _solid_perm_tensor.block(0,0,dim,dim) * _rho_GR / eta_GR;

	auto const L_pT = Eigen::MatrixXd::Zero(dim, dim);
	auto const L_px = Eigen::MatrixXd::Zero(dim, dim);

	auto const L_Tp = Eigen::MatrixXd::Zero(dim, dim);

	// TODO: add zeolite part
	auto const L_TT = Eigen::MatrixXd::Identity(dim, dim)
					  * ( _poro * lambda_F + (1.0 - _poro) * lambda_S);

	auto const L_Tx = Eigen::MatrixXd::Zero(dim, dim);

	auto const L_xp = Eigen::MatrixXd::Zero(dim, dim);
	auto const L_xT = Eigen::MatrixXd::Zero(dim, dim);

	auto const L_xx = Eigen::MatrixXd::Identity(dim, dim)
					  * _tortuosity * _poro * _rho_GR * _diffusion_coefficient_component;


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
	const double A_pp = 0.0;
	const double A_pT = 0.0;

	const double A_px = 0.0;

	const double A_Tp = 0.0;

	const double A_TT = _rho_GR * _cpG; // porosity?
	const double A_Tx = 0.0;

	const double A_xp = 0.0;
	const double A_xT = 0.0;
	const double A_xx = _rho_GR; // porosity?


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
	const double C_xx = (_poro - 1.0) * _reaction_rate;


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
	const double reaction_enthalpy = _process->getMaterials()._adsorption->get_enthalpy(_T, _p_V, _M_react);

	const double rhs_p = (_poro - 1.0) * _reaction_rate; // TODO [CL] body force term

	const double rhs_T = _rho_GR * _poro * _fluid_specific_heat_source
						 + (1.0 - _poro) * _reaction_rate * reaction_enthalpy
						 + _solid_density * (1.0 - _poro) * _solid_specific_heat_source;
						 // TODO [CL] momentum production term

	const double rhs_x = (_poro - 1.0) * _reaction_rate; // TODO [CL] what if x < 0.0


	Eigen::Vector3d rhs;
	rhs << rhs_p,
		 rhs_T,
		 rhs_x;

	return rhs;
}


void
LADataNoTpl::
preEachAssembleIntegrationPoint(const std::vector<double> &localX,
								const std::vector<double> &localSecondaryVariables,
								const VecRef &smN, const MatRef &smDNdx)
{
    auto const N = smDNdx.cols(); // number of integration points
    // auto const D = smDNdx.rows(); // global dimension

    // interpolate primary variables
    {
        double* int_pt_val[NODAL_DOF] = { &_p, &_T, &_vapour_mass_fraction };

        for (unsigned d=0; d<NODAL_DOF; ++d)
        {
            *int_pt_val[d] = 0.0;
        }

        for (unsigned d=0; d<NODAL_DOF; ++d)
        {
            for (unsigned n=0; n<N; ++n)
            {
                *int_pt_val[d] += localX[d*N+n] * smN(n);
            }
        }
    }

    // interpolate secondary variables
    {
        double* int_pt_val[NODAL_DOF_2ND] = { &_solid_density, &_reaction_rate };

        for (unsigned d=0; d<NODAL_DOF_2ND; ++d)
        {
            *int_pt_val[d] = 0.0;
        }

        for (unsigned d=0; d<NODAL_DOF_2ND; ++d)
        {
            for (unsigned n=0; n<N; ++n)
            {
                *int_pt_val[d] += localSecondaryVariables[d*N+n] * smN(n);
            }
        }
    }


#if 0
    std::cerr << "integration point values of"
                 " p=" << _p
              << " T=" << _T
              << " x=" << _vapour_mass_fraction
              << " rho_SR=" << _solid_density
              << " react_rate=" << _reaction_rate << std::endl;
#endif


    // pre-compute certain properties

    _rho_GR = fluid_density(_p, _T, _vapour_mass_fraction);
    DBUG("rho_GR = %g", _rho_GR);
    // _eta_GR = fluid_viscosity(_p, _T, _x); // used only once
    // _lambda_GR = fluid_heat_conductivity(_p, _T, _x); // used only once
    // _cpS = solid_isobaric_heat_capacity(_solid_density); // used only once
    // _H_vap = evaporation_enthalpy(_p, _T, _x); // used only once

    // _vapour_molar_fraction = Ads::Adsorption::get_molar_fraction(_vapour_mass_fraction, _M_react, _M_inert);
    _p_V = _p * Ads::Adsorption::get_molar_fraction(_vapour_mass_fraction, _M_react, _M_inert);
    DBUG("p_V = %g", _p_V);

    const double loading = Ads::Adsorption::get_loading(_solid_density, _rho_SR_dry);
    DBUG("solid_density = %g", _solid_density);


    _reaction_rate = _process->getMaterials()._adsorption->get_reaction_rate(_p_V, _T, _M_react, loading)
                     * _rho_SR_dry;
    DBUG("reaction_rate = %g", _reaction_rate);
}


void
ogs5OutMat(const LADataNoTpl::MatRef& mat)
{
    for (unsigned r=0; r<mat.rows(); ++r)
    {
        std::printf("|");

        for (unsigned c=0; c<mat.cols(); ++c)
        {
            std::printf(" %.12e", mat(r, c));
        }

        std::printf(" |\n");
    }
}

void
ogs5OutVec(const LADataNoTpl::VecRef& vec)
{
    for (unsigned r=0; r<vec.size(); ++r)
    {
        std::printf("| %.12e |\n", vec[r]);
    }
}


void
LADataNoTpl::
assembleIntegrationPoint(
		Eigen::MatrixXd* localA, Eigen::VectorXd* /*localRhs*/,
		std::vector<double> const& localX,
		std::vector<double> const& localSecondaryVariables,
		const VecRef &smN, const MatRef &smDNdx, const double smDetJ,
		const double weight)
{
    preEachAssembleIntegrationPoint(localX, localSecondaryVariables, smN, smDNdx);

    auto const N = smDNdx.cols(); // number of integration points
    auto const D = smDNdx.rows(); // global dimension: 1, 2 or 3

    assert(N*NODAL_DOF == localA->cols());
    assert(N*NODAL_DOF == _Lap->cols());

    // std::cerr << "localA:\n" << (*localA) << std::endl;
    // std::cerr << "localA block\n" << (localA->template block<2,2>(0,0)) << std::endl;
    // std::cerr << "coeff:\n" << smDNdx.transpose() * 1.0 * smDNdx * smDetJ * weight << std::endl;


    // std::cerr << "global dim:\n" << D << std::endl;

    auto const laplaceCoeffMat = getLaplaceCoeffMatrix(D);
    auto const massCoeffMat    = getMassCoeffMatrix();
    auto const advCoeffMat     = getAdvectionCoeffMatrix();
    auto const contentCoeffMat = getContentCoeffMatrix();

    Eigen::MatrixXd const velocity = Eigen::MatrixXd::Constant(D, 1, 0.0);


    DBUG("detJ = %g, weight = %g, detJ*weight = %g", smDetJ, weight, smDetJ*weight);

    std::cerr << "sm.N:\n" << smN << std::endl;
    std::cerr << "sm.dNdx:\n" << smDNdx << std::endl;

    for (unsigned r=0; r<NODAL_DOF; ++r)
    {
        for (unsigned c=0; c<NODAL_DOF; ++c)
        {
            _Lap->block(N*r, N*c, N, N).noalias() += smDetJ * weight * smDNdx.transpose() * laplaceCoeffMat.block(D*r, D*c, D, D) * smDNdx;
            _Mas->block(N*r, N*c, N, N).noalias() += smDetJ * weight * smN * massCoeffMat(r, c) * smN.transpose();
            _Adv->block(N*r, N*c, N, N).noalias() += smDetJ * weight * smN * advCoeffMat(r, c) * velocity.transpose() * smDNdx;
            _Cnt->block(N*r, N*c, N, N).noalias() += smDetJ * weight * smN * contentCoeffMat(r, c) * smN.transpose();

            // std::cerr << "vel prod: " << smDNdx.transpose() * velocity << std::endl;
#if 0
            localA->block(N*r, N*c, N, N).noalias() +=
                    (
                        smDNdx.transpose() * laplaceCoeffMat.block(D*r, D*c, D, D) * smDNdx
                        + smN * (massCoeffMat(r, c) + contentCoeffMat(r, c)) * smN.transpose()
                        + smN * advCoeffMat(r, c) * velocity.transpose() * smDNdx
                    )
                    * smDetJ * weight;
#endif

            // localA->noalias() += Lap + Mas + Adv + Cnt;

            // std::cerr << "localA block(" << N*r << "," << N*c << ")\n"
            //           << (localA->block(N*r, N*c, N, N)) << std::endl;
            // std::cerr << "coeff:\n"
            //           << smN.transpose() * massCoeffMat(r, c) * smN
            //              * smDetJ * weight << std::endl;
        }
    }

    // std::cerr << "sm.N: %s" << smN << std::endl;
    // std::cerr << "sm.dNdx: %s" << smDNdx << std::endl;
    // std::cerr << "mass coeffs:\n" << massCoeffMat << std::endl;

    auto const rhsCoeffVector = getRHSCoeffVector();

    std::printf("\n--rhs coeffs\n");
    ogs5OutVec(rhsCoeffVector);
    std::printf("\n");

    for (unsigned r=0; r<NODAL_DOF; ++r)
    {
        _rhs->block(N*r, 0, N, 1).noalias() +=
                rhsCoeffVector(r) * smN * smDetJ * weight;
    }
}


void
LADataNoTpl::preEachAssemble(const unsigned num_int_pts)
{
    if (_Lap.get() == nullptr)
    {
        _Lap.reset(new Eigen::MatrixXd(num_int_pts*NODAL_DOF, num_int_pts*NODAL_DOF));
        _Mas.reset(new Eigen::MatrixXd(num_int_pts*NODAL_DOF, num_int_pts*NODAL_DOF));
        _Adv.reset(new Eigen::MatrixXd(num_int_pts*NODAL_DOF, num_int_pts*NODAL_DOF));
        _Cnt.reset(new Eigen::MatrixXd(num_int_pts*NODAL_DOF, num_int_pts*NODAL_DOF));
        _rhs.reset(new Eigen::VectorXd(num_int_pts*NODAL_DOF));
    }


    _Lap->setZero();
    _Mas->setZero();
    _Adv->setZero();
    _Cnt->setZero();
    _rhs->setZero();
}

void
LADataNoTpl::postEachAssemble(Eigen::MatrixXd* localA, Eigen::VectorXd* localRhs,
                              Eigen::VectorXd const& oldX)
{
    localA->noalias() += *_Lap + *_Mas/_process->getMaterials()._time_step + *_Adv + *_Cnt;
    localRhs->noalias() += *_rhs
                           + *_Mas * oldX/_process->getMaterials()._time_step;

    std::printf("\nStiffness:\n");
    ogs5OutMat(*localA);
    std::printf("\n");

    std::printf("\n---Mass matrix:\n");
    ogs5OutMat(*_Mas);
    std::printf("\n");

    std::printf("---Laplacian matrix:\n");
    ogs5OutMat(*_Lap);
    std::printf("\n");

    std::printf("---Advective matrix:\n");
    ogs5OutMat(*_Adv);
    std::printf("\n");

    std::printf("---Content:\n");
    ogs5OutMat(*_Cnt);
    std::printf("\n");

    std::printf("---RHS:\n");
    ogs5OutVec(*_rhs);
    std::printf("\n");

    std::printf("---RHS2:\n");
    ogs5OutVec(*localRhs);
    std::printf("\n");
}

} // namespace TES

} // namespace ProcessLib

