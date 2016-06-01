#include "DensityHauerNaYStach.h"

namespace
{

const double c[] = {
     0.34822515198326959,     /* a0 */
     0.0078261075439284671,   /* a1 */
     0.0020337553643303471,   /* a2 */
     -1.2548427759398523e-05, /* a3 */
     -2.1060522569568419e-06, /* a4 */
     1.6607622324262219e-08,  /* a5 */
     5.3419831304985741e-10   /* a6 */
};

}

namespace Adsorption
{

double DensityHauerNaYStach::getAdsorbateDensity(const double T_Ads) const
{
    return rhoWaterHauer(T_Ads);
}


//Thermal expansivity model for water found in the works of Hauer
double DensityHauerNaYStach::getAlphaT(const double T_Ads) const
{
    // data like in python script
    const double T0 = 283.15, alpha0 = 3.781e-4; //K; 1/K

    return alpha0/(1. - alpha0 * (T_Ads-T0)); //in 1/K
}


//Characteristic curve. Return W (A)
double DensityHauerNaYStach::characteristicCurve(const double A) const
{
    double W = curvePolyfrac(c, A); //cm^3/g

    if (W < 0.0) {
        W = 0.0; // TODO [CL] debug output
    }

    return W/1.e3; //m^3/kg
}

double DensityHauerNaYStach::dCharacteristicCurve(const double A) const
{
    return dCurvePolyfrac(c, A);
}

double DensityHauerNaYStach::getEnthalpy(const double, const double, const double) const
{
    return 2055.0 * 1000.0; // J/kg value taken from Kraus' PhD thesis
}

}
