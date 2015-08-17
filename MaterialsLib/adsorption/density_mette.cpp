#include "density_mette.h"
#include "density_cook.h"

namespace Ads
{

double DensityMette::get_adsorbate_density(const double Tads) const
{
	const double T0 = 293.15;
	const double rho0 = rho_water_Dean(T0);
	const double alpha20 = alphaT_water_Dean(T0);
	return rho0 / (1. + alpha20*(Tads-T0));
}


//Thermal expansivity model for water found in the works of Hauer
double DensityMette::get_alphaT(const double Tads) const
{
	const double T0 = 293.15;
	const double alpha20 = alphaT_water_Dean(T0);
	return alpha20 / (1. + alpha20 * (Tads-T0));
}


//Characteristic curve. Return W (A)
double DensityMette::characteristic_curve(const double A) const
{
	// NaX_Mette_polyfrac_CC.pickle
	// date extracted 2015-06-23 15:38:35 file mtime 2015-06-23 15:19:26
	const double c[] = {
	    0.36340572890087813,	/* a0 */
	    -0.0013449597038375108,	/* a1 */
	    -0.0007581210111121073,	/* a2 */
	    -7.331279615575401e-08,	/* a3 */
	    5.365656973806218e-07,	/* a4 */
	    6.854673678427112e-10,	/* a5 */
	    -1.0197050219481966e-10	/* a6 */
	};

	double W = curve_polyfrac(c, A); //cm^3/g

	if (W < 0.0) {
		W = 0.0; // TODO [CL] debug output
	}

	return W/1.e3; //m^3/kg
}

}
