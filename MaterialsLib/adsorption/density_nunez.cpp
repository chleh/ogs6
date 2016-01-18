#include "density_nunez.h"

namespace
{

// NaX_Nunez_polyfrac_CC.pickle
// date extracted 2015-06-23 15:38:35 file mtime 2015-06-23 15:19:34
const double c[] = {
	0.3631900485031771,		/* a0 */
	-0.0014242280940080726,	/* a1 */
	-0.0007751726942386291,	/* a2 */
	2.1775655036811842e-08,	/* a3 */
	5.488166913667265e-07,	/* a4 */
	6.204064716725214e-10,	/* a5 */
	-1.0345385018952998e-10	/* a6 */
};

}

namespace Ads
{

double DensityNunez::get_adsorbate_density(const double Tads) const
{
	if (Tads < 273.16 or Tads > 633.15) {
		// print('Value outside admissible range for rho.');
		// return -1;
	}

	const double a[] = { 1.0644e3,-8.01905,1.445348e-2,-4.19589e-6,-4.5294e-9 };
	const double b[] = { -8.039e-3,1.8698e-5,-2.3015e-8,2.3809e-11,-1.388e-14 };
	const double u = a[0] + Tads * (a[1] + Tads * (a[2] + Tads * (a[3] + Tads * a[4]) ) );
	const double v = 1.0 + Tads * (b[0] + Tads * (b[1] + Tads * (b[2] + Tads * (b[3] + Tads * b[4]) ) ) );
	return u/v;
}


//Thermal expansivity model for water found in the works of Hauer
double DensityNunez::get_alphaT(const double Tads) const
{
	if (Tads < 273.16 or Tads > 633.15) {
		// print('Value outside admissible range for rho.');
		// return -1;
	}

	const double a[] = { 1.0644e3,-8.01905,1.445348e-2,-4.19589e-6,-4.5294e-9 };
	const double b[] = { -8.039e-3,1.8698e-5,-2.3015e-8,2.3809e-11,-1.388e-14 };
	const double u = a[0] + Tads * (a[1] + Tads * (a[2] + Tads * (a[3] + Tads * a[4]) ) );
	const double v = 1.0 + Tads * (b[0] + Tads * (b[1] + Tads * (b[2] + Tads * (b[3] + Tads * b[4]) ) ) );
	const double du = a[1] + Tads * (2.0*a[2] + Tads * (3.0*a[3] + Tads * 4.0*a[4]) );
	const double dv = b[0] + Tads * (2.0*b[1] + Tads * (3.0*b[2] + Tads * (4.0*b[3] + Tads * 5.0*b[4]) ) );
	// return - (du*v - dv*u) / u / v;
	// return (dv*u - du*v) / u / v;
	return dv/v - du/u;
}


//Characteristic curve. Return W (A)
double DensityNunez::characteristic_curve(const double A) const
{
	double W = curve_polyfrac(c, A); //cm^3/g

	if (W < 0.0) {
		W = 0.0; // TODO [CL] debug output
	}

	return W/1.e3; //m^3/kg
}

double DensityNunez::d_characteristic_curve(const double A) const
{
	return d_curve_polyfrac(c, A);
}

}
