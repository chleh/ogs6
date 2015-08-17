
#include "logog/include/logog.hpp"

#include "adsorption.h"

#include "density_legacy.h"
#include "density_100MPa.h"
#include "density_const.h"
#include "density_cook.h"
#include "density_dubinin.h"
#include "density_hauer.h"
#include "density_mette.h"
#include "density_nunez.h"


namespace Ads
{

//Saturation pressure for water used in Nunez
double Adsorption::get_equilibrium_vapour_pressure(const double Tads)
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
double Adsorption::get_evaporation_enthalpy(double Tads) //in kJ/kg
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


double Adsorption::get_molar_fraction(double xm, double M_this, double M_other)
{
	return M_other*xm/(M_other*xm + M_this*(1.0-xm));
}


double Adsorption::get_reaction_rate(const double p_Ads, const double T_Ads,
									 const double M_Ads, const double loading)
{
	const double k_rate = 6.0e-3; //to be specified

	const double A = get_potential(T_Ads, p_Ads, M_Ads);
	const double C_eq = get_adsorbate_density(T_Ads) * characteristic_curve(A);

	return k_rate * (C_eq - loading); //scaled with mass fraction
									  // this the rate in terms of loading!
}


//Evaluate adsorbtion potential A
double Adsorption::get_potential(const double T_Ads, double p_Ads, const double M_Ads) const
{
	double A = GAS_CONST * T_Ads * log(get_equilibrium_vapour_pressure(T_Ads)/p_Ads) / (M_Ads*1.e3); //in kJ/kg = J/g
	if (A < 0.0) {
		// vapour partial pressure > saturation pressure
		// A = 0.0; // TODO [CL] debug output
	}
	return A;
}


double Adsorption::get_loading(const double rho_curr, const double rho_dry)
{
	return rho_curr / rho_dry - 1.0;
}


//Calculate sorption entropy
double Adsorption::get_entropy(const double Tads, const double A) const
{
	const double epsilon = 1.0e-8;
	const double dAdlnW = 2.0*epsilon/(log(characteristic_curve(A+epsilon)) - log(characteristic_curve(A-epsilon)));
	return dAdlnW * get_alphaT(Tads);
}


//Calculate sorption enthalpy
double Adsorption::get_enthalpy(const double T_Ads, const double p_Ads, const double M_Ads) const
{
	// TODO [CL] consider using A as obtained from current loading (needs inverse CC A(W)) instead of p_Vapour, T_Vapour
	const double A = get_potential(T_Ads, p_Ads, M_Ads);

	return (get_evaporation_enthalpy(T_Ads) + A - T_Ads * get_entropy(T_Ads,A))*1000.0; //in J/kg
}


bool
stringToReactiveSystem(std::string const& name, SolidReactiveSystem& rsys_out)
{
	if (name == "Z13XBF")
		rsys_out = SolidReactiveSystem::Z13XBF;
	else if (name == "Z13XBF_100MPa")
		rsys_out = SolidReactiveSystem::Z13XBF_100MPa;
	else if (name == "Z13XBF_Const")
		rsys_out = SolidReactiveSystem::Z13XBF_Const;
	else if (name == "Z13XBF_Cook")
		rsys_out = SolidReactiveSystem::Z13XBF_Cook;
	else if (name == "Z13XBF_Dubinin")
		rsys_out = SolidReactiveSystem::Z13XBF_Dubinin;
	else if (name == "Z13XBF_Hauer")
		rsys_out = SolidReactiveSystem::Z13XBF_Hauer;
	else if (name == "Z13XBF_Mette")
		rsys_out = SolidReactiveSystem::Z13XBF_Mette;
	else if (name == "Z13XBF_Nunez")
		rsys_out = SolidReactiveSystem::Z13XBF_Nunez;
	else return false;

	return true;
}


Adsorption* Adsorption::newInstance(std::string const& rsys)
{
	SolidReactiveSystem r;
	if (stringToReactiveSystem(rsys, r)) {
		return newInstance(r);
	} else {
		ERR("unknown reactive system: %s", rsys);
		return nullptr;
	}
}


Adsorption* Adsorption::newInstance(const SolidReactiveSystem rsys)
{
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
		return nullptr;
	}
}

} // namespace Ads

