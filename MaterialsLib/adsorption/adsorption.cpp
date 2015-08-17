#include "adsorption.h"

#include "density_legacy.h"
#include "density_100MPa.h"
#include "density_const.h"
#include "density_cook.h"
#include "density_dubinin.h"
#include "density_hauer.h"
#include "density_mette.h"
#include "density_nunez.h"


const double p_min = 0.0; //minimum pressure for which we find a valid adsorption potential




namespace Ads
{

Adsorption::Adsorption()
    : rho_low(1150.0),
      M_carrier(M_N2),
      M_react(M_H2O),
      x(Eigen::VectorXd(1))
{
}


//Saturation pressure for water used in Nunez
double Adsorption::get_ps(const double Tads)
{
	//critical T and p
	const double Tc = 647.3; //K
	const double pc = 221.2e5; //Pa
	//dimensionless T
	const double Tr = Tads/Tc;
	const double theta = 1. - Tr;
	//empirical constants
	const double c[] = {-7.69123,-26.08023,-168.17065,64.23285,-118.96462,4.16717,20.97506,1.0e9,6.0};
	const double K[] = {c[0]*theta + c[1]*pow(theta,2) + c[2]*pow(theta,3) + c[3]*pow(theta,4) + c[4]*pow(theta,5),
	                    1. + c[5]*theta + c[6]*pow(theta,2)};

	const double exponent = K[0]/(K[1]*Tr) - theta/(c[7]*pow(theta,2) + c[8]);
	return pc * exp(exponent); //in Pa
}

//Evaporation enthalpy of water from Nunez
double Adsorption::get_hv(double Tads) //in kJ/kg
{
	Tads -= 273.15;
	if (Tads <= 10.){
		const double c[] = {2.50052e3,-2.1068,-3.57500e-1,1.905843e-1,-5.11041e-2,7.52511e-3,-6.14313e-4,2.59674e-5,-4.421e-7};
		double hv = 0.;
		for (size_t i=0; i< sizeof(c)/sizeof(c[0]);i++)
			hv += c[i] * pow(Tads,i);
		return hv;
	} else if (Tads <= 300.){
		const double c[] = {2.50043e3,-2.35209,1.91685e-4,-1.94824e-5,2.89539e-7,-3.51199e-9,2.06926e-11,-6.4067e-14,8.518e-17,1.558e-20,-1.122e-22};
		double hv = 0.;
		for (size_t i=0; i< sizeof(c)/sizeof(c[0]);i++)
			hv += c[i] * pow(Tads,i);
		return hv;
	} else {
		const double c[] = {2.99866e3,-3.1837e-3,-1.566964e1,-2.514e-6,2.045933e-2,1.0389e-8};
		return ((c[0] + c[2]*Tads + c[4]*pow(Tads,2))/(1. + c[1]*Tads + c[3]*pow(Tads,2) + c[5]*pow(Tads,3)));
	}
}


//evaluate specific heat capacity of adsorbate follwing Nunez
double Adsorption::get_specific_heat_capacity(const double Tads)
{
	const double c[] = {4.224,-3.716e-3,9.351e-5,-7.1786e-7,-9.1266e-9,2.69247e-10,-2.773104e-12,1.553177e-14,-4.982795e-17,8.578e-20,-6.12423e-23};
	double cp = 0.;
	for (unsigned i=0; i< sizeof(c)/sizeof(c[0]);i++)
		cp += c[i] * pow(Tads,i);
	return cp; //kJ/(kg*K)
}



double Adsorption::get_mole_fraction(double xm) const
{
	return M_carrier*xm/(M_carrier*xm + M_react*(1.0-xm));
}


void Adsorption::calculate_qR()
{
	//Convert mass fraction into mole fraction
	const double mol_frac_react = get_mole_fraction(x_react);

	//partial pressure
	p_r_g = std::max(mol_frac_react * p_gas * 1.0e5, p_min); //avoid illdefined log, gas pressure in Pa
	set_sorption_equilibrium();
	const double dCdt = Z13XBF_adsorption();
	qR = rho_low * dCdt;
}


void Adsorption::set_rho_s(double new_rho_s)
{
	rho_s = new_rho_s;
}


double Adsorption::get_qR() const
{
	return qR;
}


void Adsorption::get_x(Eigen::VectorXd& output_x) const
{
	output_x = x;
}


double Adsorption::Z13XBF_adsorption() const
{
	const double k_rate = 6.0e-3; //to be specified
	const double C = rho_s/rho_low - 1.; //current degree of loading
	const double dCdt = k_rate * (C_eq - C); //scaled with mass fraction
	if (dCdt > 0. && p_r_g < p_min*1.01)
		return 0.;
	//automatic time stepping should be used instead of the following.
	//else if (dCdt > 0.) {
	//	double dens_guess = p_r_g*COMP_MOL_MASS_WATER/(R*T); //vapor density guess
	//	//const double max_rate = dens_guess/(rho_low*dt);
	//	//dCdt = min(dCdt,max_rate);
	//}
	return dCdt;
}


//determine equilibrium loading according to Dubinin
void Adsorption::set_sorption_equilibrium()
{
	//determine adsorption potential
	const double A = get_potential(T_s, p_r_g);
	//determine adsorbed volume
	const double W = characteristic_curve(A);
	//determine equilibrium loading
	// TODO [CL] eliminate this member variable
	C_eq = W * get_adsorbate_density(T_s); //kg/kg
}


//Evaluate adsorbtion potential A
double Adsorption::get_potential(const double Tads, double pads) const
{
	pads = std::max(pads, p_min);
	double A = GAS_CONST * Tads * log(get_ps(Tads)/pads) / (M_react*1.e3); //in kJ/kg = J/g
	if (A < 0.0) {
		// vapour partial pressure > saturation pressure
		// A = 0.0; // TODO [CL] debug output
	}
	return A;
}


void Adsorption::update_param(double T_solid,
                              double p_gas,
                              double x_reactive,
                              double rho_s_initial)
{
	this->T_s       = T_solid;
	// TODO: combine to partial pressure of reactive component
	this->p_gas     = p_gas; // should be in unit bar
	this->x_react   = x_reactive;

	this->rho_s     = rho_s_initial;
	this->x(0)      = rho_s_initial;
}


void Adsorption::eval(double /*t*/, Eigen::VectorXd const& y, Eigen::VectorXd &dydx)
{
	assert( y.size() == dydx.size() );

	set_rho_s( y(0) );
	calculate_qR();
	dydx(0) = get_qR();

}


//Calculate sorption entropy
double Adsorption::get_entropy(const double Tads, const double A) const
{
	const double epsilon = 1.0e-8;
	const double dAdlnW = 2.0*epsilon/(log(characteristic_curve(A+epsilon)) - log(characteristic_curve(A-epsilon)));
	return dAdlnW * get_alphaT(Tads);
}


//Calculate sorption enthalpy
double Adsorption::get_enthalpy(const double Tads, const double pads) const
{
	// TODO [CL] consider using A as obtained from current loading (needs inverse CC A(W)) instead of p_Vapour, T_Vapour
	const double A = get_potential(Tads, pads);

	return (get_hv(Tads) + A - Tads * get_entropy(Tads,A))*1000.0; //in J/kg
}


Adsorption* Adsorption::newInstance(const SolidReactiveSystem rsys)
{
	// using ::FiniteElement;

	switch (rsys)
	{
	case SolidReactiveSystem::Z13XBF:
		return new DensityLegacy();
	case SolidReactiveSystem::Z13XBF_100MPa:
		return new Density100MPa();
	case SolidReactiveSystem::Z13XBF_Const:
		return new DensityConst();
	case SolidReactiveSystem::Z13XBF_Cook:
		return new DensityCook();
	case SolidReactiveSystem::Z13XBF_Dubinin:
		return new DensityDubinin();
	case SolidReactiveSystem::Z13XBF_Hauer:
		return new DensityHauer();
	case SolidReactiveSystem::Z13XBF_Mette:
		return new DensityMette();
	case SolidReactiveSystem::Z13XBF_Nunez:
		return new DensityNunez();
	default:
		return NULL;
	}
}

}

