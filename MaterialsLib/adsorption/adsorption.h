#ifndef ADSORPTION_H
#define ADSORPTION_H

#include "Eigen/Eigen"

namespace Ads
{


const double GAS_CONST = 8.3144621;

const double M_N2  = 0.028;
const double M_H2O = 0.018;

enum class SolidReactiveSystem
{
	Z13XBF,
	Z13XBF_Const,
	Z13XBF_Hauer,
	Z13XBF_Mette,
	Z13XBF_Nunez,
	Z13XBF_Cook,
	Z13XBF_Dubinin,
	Z13XBF_100MPa,
	Inert
};

class Adsorption
{
public:
// static members
	static Adsorption* newInstance(SolidReactiveSystem rsys);
	static Adsorption* newInstance(std::string const& rsys);

	// TODO [CL] move those four methods to water properties class
	static double get_evaporation_enthalpy(const double Tads);
	static double get_equilibrium_vapour_pressure(const double Tads);
	static double get_specific_heat_capacity(const double Tads); // TODO [CL] why unused?
	static double get_molar_fraction(double xm, double M_this, double M_other);

	static double get_loading(const double rho_curr, const double rho_dry);

// virtual members:
	virtual ~Adsorption() {}

	virtual double get_enthalpy(const double T_Ads, const double p_Ads, const double M_Ads) const;
	virtual double get_reaction_rate(const double p_Ads, const double T_ads,
									 const double M_Ads, const double loading) const;

protected:
	virtual double get_adsorbate_density(const double Tads) const = 0;
	virtual double get_alphaT(const double Tads) const = 0;
	virtual double characteristic_curve(const double A) const = 0;

private:
// non-virtual members
	double get_potential(const double T_Ads, const double p_Ads, const double M_Ads) const;
	double get_entropy(const double Tads, const double A) const;
};


inline double curve_polyfrac(const double* coeffs, const double x)
{
	return ( coeffs[0] + coeffs[2] * x + coeffs[4] * pow(x,2) + coeffs[6] * pow(x,3) )
	        / ( 1.0 + coeffs[1] * x + coeffs[3] * pow(x,2) + coeffs[5] * pow(x,3) );

	// Apparently, even pow and std::pow are different
	// return ( coeffs[0] + coeffs[2] * x + coeffs[4] * std::pow(x,2) + coeffs[6] * std::pow(x,3) )
	//         / ( 1.0 + coeffs[1] * x + coeffs[3] * std::pow(x,2) + coeffs[5] * std::pow(x,3) );

	// Analytically the same, but numerically quite different
	// return ( coeffs[0] + x * ( coeffs[2] + x * (coeffs[4] + x * coeffs[6] ) ) )
	//         / ( 1.0 + x * ( coeffs[1] + x * (coeffs[3] + x * coeffs[5] ) ) );

	// Analytically the same, but numerically quite different
	// return ( coeffs[0] + x * coeffs[2] + x*x * coeffs[4] + x*x*x * coeffs[6] )
	//        / ( 1.0 + x * coeffs[1] + x*x * coeffs[3] + x*x*x * coeffs[5] );
}

}

#endif // ADSORPTION_H
