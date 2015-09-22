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

#include <typeinfo>

#include "logog/include/logog.hpp"

#include "TESFEM-notpl.h"

#include "MathLib/Nonlinear/Root1D.h"
#include "NumLib/Function/Interpolation.h"


const double GAS_CONST = 8.3144621;

enum class MatOutType { OGS5, PYTHON };

const MatOutType MATRIX_OUTPUT_FORMAT = MatOutType::PYTHON;



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
	// const double R = GAS_CONST;
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
getMassCoeffMatrix(const unsigned int_pt)
{
	double dxn_dxm = _AP->_M_inert * _AP->_M_react
					 / square(_AP->_M_inert * _vapour_mass_fraction
							  + _AP->_M_react * (1.0 - _vapour_mass_fraction));

	const double M_pp = _AP->_poro/_p * _rho_GR;
	const double M_pT = -_AP->_poro/_T *  _rho_GR;
	const double M_px = (_AP->_M_react-_AP->_M_inert) * _p
						/ (GAS_CONST * _T) * dxn_dxm * _AP->_poro;

	const double M_Tp = -_AP->_poro;
	const double M_TT = _AP->_poro * _rho_GR * _AP->_cpG // TODO: vapour heat capacity
						+ (1.0-_AP->_poro) * _solid_density[int_pt] * _AP->_cpS; // TODO: adsorbate heat capacity
	const double M_Tx = 0.0;

	const double M_xp = 0.0;
	const double M_xT = 0.0;
	const double M_xx = _AP->_poro * _rho_GR;


	Eigen::Matrix3d M;
	M << M_pp, M_pT, M_px,
		 M_Tp, M_TT, M_Tx,
		 M_xp, M_xT, M_xx;

	return M;
}


Eigen::MatrixXd
LADataNoTpl::
getLaplaceCoeffMatrix(const unsigned /*int_pt*/, const unsigned dim)
{
	const double eta_GR = fluid_viscosity(_p, _T, _vapour_mass_fraction);

	const double lambda_F = fluid_heat_conductivity(_p, _T, _vapour_mass_fraction);
	const double lambda_S = _AP->_solid_heat_cond;

	// TODO: k_rel
	Eigen::MatrixXd L_pp = _AP->_solid_perm_tensor.block(0,0,dim,dim) * _rho_GR / eta_GR;

	Eigen::MatrixXd L_pT = Eigen::MatrixXd::Zero(dim, dim);
	Eigen::MatrixXd L_px = Eigen::MatrixXd::Zero(dim, dim);

	Eigen::MatrixXd L_Tp = Eigen::MatrixXd::Zero(dim, dim);

	// TODO: add zeolite part
	Eigen::MatrixXd L_TT = Eigen::MatrixXd::Identity(dim, dim)
					  * ( _AP->_poro * lambda_F + (1.0 - _AP->_poro) * lambda_S);

	Eigen::MatrixXd L_Tx = Eigen::MatrixXd::Zero(dim, dim);

	Eigen::MatrixXd L_xp = Eigen::MatrixXd::Zero(dim, dim);
	Eigen::MatrixXd L_xT = Eigen::MatrixXd::Zero(dim, dim);

	Eigen::MatrixXd L_xx = Eigen::MatrixXd::Identity(dim, dim)
						   * (_AP->_tortuosity * _AP->_poro * _rho_GR
							  * _AP->_diffusion_coefficient_component);

	Eigen::MatrixXd L(dim*3, dim*3);

	L.block(    0,     0, dim, dim) = L_pp;
	L.block(    0,   dim, dim, dim) = L_pT;
	L.block(    0, 2*dim, dim, dim) = L_px;

	L.block(  dim,     0, dim, dim) = L_Tp;
	L.block(  dim,   dim, dim, dim) = L_TT;
	L.block(  dim, 2*dim, dim, dim) = L_Tx;

	L.block(2*dim,     0, dim, dim) = L_xp;
	L.block(2*dim,   dim, dim, dim) = L_xT;
	L.block(2*dim, 2*dim, dim, dim) = L_xx;

	return L;
}


Eigen::Matrix3d
LADataNoTpl::
getAdvectionCoeffMatrix(const unsigned /*int_pt*/)
{
	const double A_pp = 0.0;
	const double A_pT = 0.0;

	const double A_px = 0.0;

	const double A_Tp = 0.0;

	const double A_TT = _rho_GR * _AP->_cpG; // porosity?
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
getContentCoeffMatrix(const unsigned /*int_pt*/)
{
	const double C_pp = 0.0;
	const double C_pT = 0.0;

	const double C_px = 0.0;

	const double C_Tp = 0.0;

	const double C_TT = 0.0;
	const double C_Tx = 0.0;

	const double C_xp = 0.0;
	const double C_xT = 0.0;
	const double C_xx = (_AP->_poro - 1.0) * _qR;


	Eigen::Matrix3d C;
	C << C_pp, C_pT, C_px,
		 C_Tp, C_TT, C_Tx,
		 C_xp, C_xT, C_xx;

	return C;
}


Eigen::Vector3d
LADataNoTpl::
getRHSCoeffVector(const unsigned int_pt)
{
	const double reaction_enthalpy = _AP->_adsorption->get_enthalpy(_p_V, _T, _AP->_M_react);

	const double rhs_p = (_AP->_poro - 1.0) * _qR; // TODO [CL] body force term

	const double rhs_T = _rho_GR * _AP->_poro * _AP->_fluid_specific_heat_source
						 + (1.0 - _AP->_poro) * _qR * reaction_enthalpy
						 + _solid_density[int_pt] * (1.0 - _AP->_poro) * _AP->_solid_specific_heat_source;
						 // TODO [CL] momentum production term

	const double rhs_x = (_AP->_poro - 1.0) * _qR; // TODO [CL] what if x < 0.0


	Eigen::Vector3d rhs;
	rhs << rhs_p,
		 rhs_T,
		 rhs_x;

	return rhs;
}


void LADataNoTpl::
initNewTimestep(const unsigned int_pt, const std::vector<double> &/*localX*/)
{
    if (_AP->_iteration_in_current_timestep == 0)
    {
        // first get reaction rate from the given kinetics
        // this does not consider that vapour is actually sucked up into the zeolite
        const double loading = Ads::Adsorption::get_loading(_solid_density_prev_ts[int_pt], _AP->_rho_SR_dry);
        auto const react_rate_R = _AP->_adsorption->get_reaction_rate(_p_V, _T, _AP->_M_react, loading)
                                  * _AP->_rho_SR_dry;

        // calculate density change
        const double delta_rhoS = react_rate_R * _AP->_delta_t * (1.0 - _AP->_poro);
        const double delta_rhoV   = - delta_rhoS;
        const double rho_V = _AP->_M_react * _p_V / GAS_CONST / _T * _AP->_poro;

        if (-delta_rhoV > rho_V)
        {
            // in this case with the model only considering adsorption kinetics, the zeolite will adsorb more water than
            // there actually is ==> limit adsorption, use equilibrium reaction

            // function describing local equilibrium between vapour and zeolite loading
            // temperature is assumed to be constant
            auto f = [this, loading](double pV) -> double
            {
                // pV0 := _pV
                const double C_eq = _AP->_adsorption->get_equilibrium_loading(pV, _T, _AP->_M_react);
                return (pV - _p_V) * _AP->_M_react / GAS_CONST / _T * _AP->_poro
                        + (1.0-_AP->_poro) * (C_eq - loading) * _AP->_rho_SR_dry;
            };

            // range where to search for roots of f
            const double C_eq0 = _AP->_adsorption->get_equilibrium_loading(_p_V, _T, _AP->_M_react);
            const double limit = (C_eq0 > loading) ? 1e-8 : _p; // TODO [CL] upper limit to equilibrium vapour pressure

            // search for roots
            auto rf = MathLib::Nonlinear::makeRegulaFalsi<MathLib::Nonlinear::Pegasus>(f, _p_V, limit);
            rf.step(3);

            // set vapour pressure
            const double pV = rf.get_result();
            const double delta_pV = pV - _p_V;
            _p += delta_pV;
            _p_V = pV;

            // set solid density
            const double delta_rhoV = delta_pV * _AP->_M_react / GAS_CONST / _T * _AP->_poro;
            const double delta_rhoSR = delta_rhoV / (_AP->_poro - 1.0);
            _reaction_rate[int_pt] = delta_rhoSR / _AP->_delta_t;
            _solid_density[int_pt] = _solid_density_prev_ts[int_pt] + delta_rhoSR;

            DBUG("pV: %14.7g, delta pV: %14.7g, rhoSR: %14.7g, delta_rhoSR: %14.7g", _p_V, delta_pV, _solid_density[int_pt], delta_rhoSR);
        }
        else
        {
            // in this case considering only the adsorption kintetics is sufficient,
            // and may be also more correct, since it really describes a kinetic.

            _reaction_rate[int_pt] = react_rate_R;
            _solid_density[int_pt] = _solid_density_prev_ts[int_pt] + react_rate_R * _AP->_delta_t;
        }

        _qR = _reaction_rate[int_pt];
    }

    return;


    const double loading = Ads::Adsorption::get_loading(_solid_density[int_pt], _AP->_rho_SR_dry);

    auto const dCdt = _AP->_adsorption->get_reaction_rate(_p_V, _T, _AP->_M_react, loading);

    auto              qR = dCdt * _AP->_rho_SR_dry;
    auto const delta_rho = qR * _AP->_delta_t;          // solid density change in this timestep
    auto const     rho_V = M_H2O * _p_V / GAS_CONST / _T; // vapour density

    const double rate_factor = 1.0;

    if (delta_rho > rate_factor * rho_V) {
        // all vapour will be sucked up in this time step
        // => limit reaction rate
        // DBUG("limiting reaction rate %g --> %g", qR, qR * rho_V / delta_rho);
        qR *= rate_factor * rho_V / delta_rho;
    }

    if (qR > 0.0 && _p_V < 100)
    {
        // qR = 0.0;
    }

    /*
    // average reaction rate between this and the previous iteration
    double average_qR = (_AP->_iteration_in_current_timestep < 2)
            ? qR
            : (0.5 * (qR + _reaction_rate[int_pt]));

    _qR = average_qR; // this now only works if reaction rate is calculated in every timestep
    */

    _qR = qR;

    /*
    if ((average_qR < 0.0) != (_reaction_rate[int_pt] < 0.0)) {
        // old and new reaction rates have different sign.
        // let the reaction take a break of one iteration
        average_qR = 0.0;
    }
    */

    // dCdt = qR / _AP->_rho_SR_dry;

    _reaction_rate[int_pt] = qR;
    _solid_density[int_pt] = _solid_density_prev_ts[int_pt] + _qR * _AP->_delta_t;
    // _solid_density[int_pt] = _AP->_rho_SR_dry * (1.0 + C_next);
}


void
LADataNoTpl::
preEachAssembleIntegrationPoint(
        const unsigned int_pt,
        const std::vector<double> &localX,
        const VecRef &smN, const MatRef& /*smDNdx*/)
{
    std::array<double*, NODAL_DOF> int_pt_val = { &_p, &_T, &_vapour_mass_fraction };

    NumLib::shapeFunctionInterpolate(localX, smN, int_pt_val);

    //*
    if (_p < 1.0) _p = 1.0;
    if (_T < 274.0) _T = 274.0;
    else if (_T > 600.0) _T = 600.0;
    if (_vapour_mass_fraction < 1e-6) _vapour_mass_fraction = 1e-6;
    else if (_vapour_mass_fraction > 1.0 - 1e-6) _vapour_mass_fraction = 1.0 - 1e-6;
    //*/

    assert(_p > 0.0);
    assert(_T > 0.0);
    assert(0.0 <= _vapour_mass_fraction && _vapour_mass_fraction <= 1.0);

    // pre-compute certain properties
    _rho_GR = fluid_density(_p, _T, _vapour_mass_fraction);
    _p_V = _p * Ads::Adsorption::get_molar_fraction(_vapour_mass_fraction, _AP->_M_react, _AP->_M_inert);

    // if (_p_V <= 0.0) _p_V = std::numeric_limits<double>::epsilon();

    if (true || _AP->_iteration_in_current_timestep == 0) {
        initNewTimestep(int_pt, localX);
    }
}


void
ogs5OutMat(const LADataNoTpl::MatRef& mat)
{
    for (unsigned r=0; r<mat.rows(); ++r)
    {
        switch (MATRIX_OUTPUT_FORMAT)
        {
        case MatOutType::OGS5:
            if (r!=0) std::printf("\n");
            std::printf("|");
            break;
        case MatOutType::PYTHON:
            if (r!=0) std::printf(",\n");
            std::printf("[");
            break;
        }

        for (unsigned c=0; c<mat.cols(); ++c)
        {
            switch (MATRIX_OUTPUT_FORMAT)
            {
            case MatOutType::OGS5:
                std::printf(" %.16e", mat(r, c));
                break;
            case MatOutType::PYTHON:
                if (c!=0) std::printf(",");
                std::printf(" %23.16g", mat(r, c));
                break;
            }

        }

        switch (MATRIX_OUTPUT_FORMAT)
        {
        case MatOutType::OGS5:
            std::printf(" | ");
            break;
        case MatOutType::PYTHON:
            std::printf(" ]");
            break;
        }
    }
    std::printf("\n");
}

void
ogs5OutVec(const LADataNoTpl::VecRef& vec)
{
    for (unsigned r=0; r<vec.size(); ++r)
    {
        switch (MATRIX_OUTPUT_FORMAT)
        {
        case MatOutType::OGS5:
            if (r!=0) std::printf("\n");
            std::printf("| %.16e | ", vec[r]);
            break;
        case MatOutType::PYTHON:
            if (r!=0) std::printf(",\n");
            std::printf("[ %23.16g ]", vec[r]);
            break;
        }
    }
    std::printf("\n");
}


std::vector<double> const&
LADataNoTpl::
getIntegrationPointValues(SecondaryVariables var) const
{
    switch (var)
    {
    case SecondaryVariables::REACTION_RATE:
        return _reaction_rate;
        break;
    case SecondaryVariables::SOLID_DENSITY:
        return _solid_density;
        break;
    case SecondaryVariables::VELOCITY_X:
        return _velocity[0];
    case SecondaryVariables::VELOCITY_Y:
        assert(_velocity.size() >= 2);
        return _velocity[1];
    case SecondaryVariables::VELOCITY_Z:
        assert(_velocity.size() >= 3);
        return _velocity[2];
    }

    // TODO: error!
    return _reaction_rate;
}


void
LADataNoTpl::
assembleIntegrationPoint(unsigned integration_point,
                         Eigen::MatrixXd* localA, Eigen::VectorXd* /*localRhs*/,
                         std::vector<double> const& localX,
                         const VecRef &smN, const MatRef &smDNdx, const double smDetJ,
                         const double weight)
{
    preEachAssembleIntegrationPoint(integration_point, localX, smN, smDNdx);

    auto const N = smDNdx.cols(); // number of integration points
    auto const D = smDNdx.rows(); // global dimension: 1, 2 or 3

    assert(N*NODAL_DOF == localA->cols());
    assert(N*NODAL_DOF == _Lap->cols());

    auto const laplaceCoeffMat = getLaplaceCoeffMatrix(integration_point, D);
    assert(laplaceCoeffMat.cols() == D*NODAL_DOF);
    auto const massCoeffMat    = getMassCoeffMatrix(integration_point);
    auto const advCoeffMat     = getAdvectionCoeffMatrix(integration_point);
    auto const contentCoeffMat = getContentCoeffMatrix(integration_point);


    // calculate velocity
    assert((unsigned) smDNdx.rows() == _velocity.size() && (unsigned) smDNdx.cols() == _velocity[0].size());

    // using auto for the type went terribly wrong!
    // calculating grad_p not separately also went wrong!
    Eigen::VectorXd const grad_p = smDNdx * Eigen::Map<const Eigen::VectorXd>(localX.data(), N);
    assert(grad_p.size() == D);
    Eigen::VectorXd const velocity = - laplaceCoeffMat.block(0, 0, D, D) * grad_p
                                     / _rho_GR;
    assert(velocity.size() == D);

    for (unsigned d=0; d<D; ++d)
    {
        _velocity[d][integration_point] = velocity[d];
    }

    Eigen::VectorXd const detJ_w_N = smDetJ * weight * smN;
    Eigen::MatrixXd const detJ_w_N_NT = detJ_w_N * smN.transpose();
    assert(detJ_w_N_NT.rows() == N && detJ_w_N_NT.cols() == N);

    Eigen::MatrixXd const vT_dNdx = velocity.transpose() * smDNdx;
    assert(vT_dNdx.cols() == N && vT_dNdx.rows() == 1);
    Eigen::MatrixXd const detJ_w_N_vT_dNdx = detJ_w_N * vT_dNdx;
    assert(detJ_w_N_vT_dNdx.rows() == N && detJ_w_N_vT_dNdx.cols() == N);

    for (unsigned r=0; r<NODAL_DOF; ++r)
    {
        for (unsigned c=0; c<NODAL_DOF; ++c)
        {
            Eigen::MatrixXd tmp = smDetJ * weight * smDNdx.transpose();
            assert(tmp.cols() == D && tmp.rows() == N);
            tmp *= laplaceCoeffMat.block(D*r, D*c, D, D);
            assert(tmp.cols() == D && tmp.rows() == N);
            tmp *= smDNdx;
            assert(tmp.cols() == N && tmp.rows() == N);

            _Lap->block(N*r, N*c, N, N).noalias() += tmp;
            _Mas->block(N*r, N*c, N, N).noalias() += detJ_w_N_NT      * massCoeffMat(r, c);
            _Adv->block(N*r, N*c, N, N).noalias() += detJ_w_N_vT_dNdx * advCoeffMat(r, c);
            _Cnt->block(N*r, N*c, N, N).noalias() += detJ_w_N_NT      * contentCoeffMat(r, c);
        }
    }

    auto const rhsCoeffVector = getRHSCoeffVector(integration_point);

    for (unsigned r=0; r<NODAL_DOF; ++r)
    {
        _rhs->block(N*r, 0, N, 1).noalias() +=
                rhsCoeffVector(r) * smN * smDetJ * weight;
    }
}


void
LADataNoTpl::init(const unsigned num_int_pts, const unsigned dimension)
{
    _solid_density.resize(num_int_pts, _AP->_initial_solid_density);
    _solid_density_prev_ts.resize(num_int_pts, _AP->_initial_solid_density);

    _reaction_rate.resize(num_int_pts);
    _reaction_rate_prev_ts.resize(num_int_pts);

    // _velocity.resize(num_int_pts, dimension);
    _velocity.resize(dimension);
    for (auto& v : _velocity) v.resize(num_int_pts);

    _Lap.reset(new Eigen::MatrixXd(num_int_pts*NODAL_DOF, num_int_pts*NODAL_DOF));
    _Mas.reset(new Eigen::MatrixXd(num_int_pts*NODAL_DOF, num_int_pts*NODAL_DOF));
    _Adv.reset(new Eigen::MatrixXd(num_int_pts*NODAL_DOF, num_int_pts*NODAL_DOF));
    _Cnt.reset(new Eigen::MatrixXd(num_int_pts*NODAL_DOF, num_int_pts*NODAL_DOF));
    _rhs.reset(new Eigen::VectorXd(num_int_pts*NODAL_DOF));

    _Lap->setZero();
    _Mas->setZero();
    _Adv->setZero();
    _Cnt->setZero();
    _rhs->setZero();
}


void
LADataNoTpl::preEachAssemble()
{
    if (_AP->_iteration_in_current_timestep == 0)
    {
        _solid_density_prev_ts = _solid_density;
        _reaction_rate_prev_ts = _reaction_rate;
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
    localA->noalias() += *_Lap + *_Mas/_AP->_delta_t + *_Adv + *_Cnt;
    localRhs->noalias() += *_rhs
                           + *_Mas * oldX/_AP->_delta_t;

    if (_AP->_output_element_matrices)
    {
        std::puts("### Element: ?");

        std::puts("---Velocity of water");
        for (auto const& vs : _velocity)
        {
            std::printf("| ");
            for (auto v : vs)
            {
                std::printf("%23.16e ", v);
            }
            std::printf("|\n");
        }

        std::printf("\nStiffness: \n");
        ogs5OutMat(*localA);
        std::printf("\n");

        std::printf("\n---Mass matrix: \n");
        ogs5OutMat(*_Mas);
        std::printf("\n");

        std::printf("---Laplacian matrix: \n");
        ogs5OutMat(*_Lap);
        std::printf("\n");

        std::printf("---Advective matrix: \n");
        ogs5OutMat(*_Adv);
        std::printf("\n");

        std::printf("---Content: \n");
        ogs5OutMat(*_Cnt);
        std::printf("\n");

        std::printf("---RHS: \n");
        ogs5OutVec(*localRhs);
        std::printf("\n");
    }
}

} // namespace TES

} // namespace ProcessLib

