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
	Z13XBF_100MPa
};

class Adsorption
{
// static members
public:
	static Adsorption* newInstance(SolidReactiveSystem rsys);

	static double get_hv(const double Tads);
	static double get_ps(const double Tads);
	static double get_specific_heat_capacity(const double Tads); // TODO [CL] why unused?


// virtual members:
public:
	virtual ~Adsorption() {}

	virtual double get_adsorbate_density(const double Tads) const = 0;
	virtual double get_alphaT(const double Tads) const = 0;
	virtual double characteristic_curve(const double A) const = 0;


// "normal" members
public:
	double get_mole_fraction(double xm) const;
	double get_enthalpy(const double Tads, const double pads) const;
	double get_potential(const double Tads, double pads) const;
	double get_entropy(const double Tads, const double A) const;

	Adsorption();

	void update_param(double T_solid,
	                  double p_gas,
	                  double w_water,
	                  double rho_s_initial);

	void eval(double t, Eigen::VectorXd const& y, Eigen::VectorXd &dydx);

private:
	double get_qR() const;
	void get_x(Eigen::VectorXd& output_x) const;
	double Z13XBF_adsorption() const; // TODO [CL] rename

	void calculate_qR();
	void set_rho_s(double new_rho_s);
	void set_sorption_equilibrium();

	// TODO [CL] reduce number of variables
	double rho_s;    // solid phase density
	double p_gas;    // [bar] gas phase pressure;
	double p_r_g;    // [Pa]  pressure of H2O on gas phase;
	double T_s;      // solid phase temperature;
	double qR;       // rate of solid density change;
	double x_react;    // mass fraction of water in gas phase;
	const double rho_low; //lower density limit
	const double M_carrier; //inert component molar mass
	const double M_react;   //reactive component molar mass
	double C_eq; //equilibrium loading of sorbens

	Eigen::VectorXd x;
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
