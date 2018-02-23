/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "DensityHauerNaYStach.h"
#include "Adsorption.h"
#include "DensityHauer.h"
#include "MaterialLib/PhysicalConstant.h"

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
const double DensityHauerNaYStach::M_Ads =
    MaterialLib::PhysicalConstant::MolarMass::Water;

double DensityHauerNaYStach::getAdsorbateDensity(const double T_Ads)
{
    return DensityHauer::getAdsorbateDensity(T_Ads);
}

// Thermal expansivity model for water found in the works of Hauer
double DensityHauerNaYStach::getAlphaT(const double T_Ads)
{
    return DensityHauer::getAlphaT(T_Ads);
}

// Characteristic curve. Return W (A)
double DensityHauerNaYStach::characteristicCurve(const double A)
{
    double W = curvePolyfrac(c, A);  // cm^3/g

    if (W < 0.0)
    {
        W = 0.0;  // TODO [CL] debug output
    }

    return W / 1.e3;  // m^3/kg
}

double DensityHauerNaYStach::dCharacteristicCurve(const double A)
{
    return dCurvePolyfrac(c, A) / 1000.0;
}

double DensityHauerNaYStach::getEnthalpy(double, double)
{
    return 2055.0 * 1000.0;  // J/kg value taken from Kraus' PhD thesis
}
}
