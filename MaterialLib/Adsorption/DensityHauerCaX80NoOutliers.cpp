/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "DensityHauerCaX80NoOutliers.h"
#include "Adsorption.h"
#include "MaterialLib/PhysicalConstant.h"

namespace
{
// CaX(80) no desorption, no outliers
// out-appl-en-conf-2016/TN060_ohne_desorp_ohne_null_Hauer_polyfrac_CC.pickle
// date extracted 2016-06-27 10:03:09 file mtime 2016-06-27 09:43:46
const double c[] = {
    0.3345380947992251,     /* a0 */
    -0.0011580333609898259, /* a1 */
    -0.0005253928557744893, /* a2 */
    4.38211443758908e-07,   /* a3 */
    2.933844269139041e-07,  /* a4 */
    8.052413297861738e-11,  /* a5 */
    -5.443882255079396e-11  /* a6 */
};
}

namespace Adsorption
{
const double DensityHauerCaX80NoOutliers::M_Ads =
    MaterialLib::PhysicalConstant::MolarMass::Water;

double DensityHauerCaX80NoOutliers::getAdsorbateDensity(const double T_Ads)
{
    return DensityHauer::getAdsorbateDensity(T_Ads);
}

double DensityHauerCaX80NoOutliers::getAlphaT(const double T_Ads)
{
    return DensityHauer::getAlphaT(T_Ads);
}

// Characteristic curve. Return W (A)
double DensityHauerCaX80NoOutliers::characteristicCurve(const double A)
{
    double W = curvePolyfrac(c, A);  // cm^3/g

    if (W < 0.0)
    {
        W = 0.0;  // TODO [CL] debug output
    }

    return W / 1.e3;  // m^3/kg
}

double DensityHauerCaX80NoOutliers::dCharacteristicCurve(const double A)
{
    return dCurvePolyfrac(c, A) / 1000.0;
}
}
