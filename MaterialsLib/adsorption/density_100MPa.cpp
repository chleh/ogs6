#include "density_100MPa.h"

namespace Ads
{

double Density100MPa::get_adsorbate_density(const double Tads) const
{
	return -0.0013*Tads*Tads + 0.3529*Tads + 1049.2;
}


//Thermal expansivity model for water found in the works of Hauer
double Density100MPa::get_alphaT(const double Tads) const
{
	const double rho    = -0.0013*Tads*Tads+0.3529*Tads+1049.2;
	const double drhodT = -0.0026*Tads + 0.3529;

	return - drhodT / rho;
}


//Characteristic curve. Return W (A)
double Density100MPa::characteristic_curve(const double A) const
{
	// NaX_HighP_polyfrac_CC.pickle
	// date extracted 2015-06-23 15:38:35 file mtime 2015-06-23 15:19:57
	const double c[] = {
	    0.3490302932983226,		/* a0 */
	    -0.0014061345691831226,	/* a1 */
	    -0.0007399303393402753,	/* a2 */
	    5.129318840267485e-09,	/* a3 */
	    5.243619689772646e-07,	/* a4 */
	    6.347011955956523e-10,	/* a5 */
	    -9.919599580166727e-11	/* a6 */
	};

	double W = curve_polyfrac(c, A); //cm^3/g

	if (W < 0.0) {
		W = 0.0; // TODO [CL] debug output
	}

	return W/1.e3; //m^3/kg
}

}

