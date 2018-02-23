/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "DensityHauer15CaCl2CaX.h"

#include "MaterialLib/PhysicalConstant.h"

#include "Adsorption.h"
#include "DensityHauer.h"

namespace
{
#if 0
// 15CaCl2/CaX ion exchanged Temperatures binned, no outliers
// $ ./get-all-cc-params.py CaCl2_CaX_Kerskes/cax+15cc__Hauer_polyfrac_CC.pickle
// CaCl2_CaX_Kerskes/cax+15cc__Hauer_polyfrac_CC.pickle
// date extracted 2017-04-13 11:50:09 file mtime 2017-04-13 11:46:51
const double c[] = {
    1.0173792852763517,      /* a0 */
    0.017506926303366886,    /* a1 */
    -0.002087006695623531,   /* a2 */
    -2.8097303437857505e-05, /* a3 */
    7.271150618497568e-06,   /* a4 */
    4.115906544344751e-08,   /* a5 */
    -2.4576306320368285e-09  /* a6 */
};
#endif

// 15CaCl2/CaX ion exchanged Temperatures binned, no outliers
// Resiudal weighted:
//
// e = np.ones_like(x) # fake error value
// idcs = np.where(x<600.0)
// print(idcs)
// np.put(e, idcs, 0.1 + 0.9 * (x[idcs]/600.0)**15)
//
// return (model - data) / e
//
// PFWeighted_cax+15cc_Hauer_polyfrac_CC.pickle
// date extracted 2017-05-17 11:53:49 file mtime 2017-05-17 11:47:22
const double c[] = {
    0.9487630041128019,     /* a0 */
    0.012093168940018295,   /* a1 */
    -0.0033190906417513147, /* a2 */
    -1.839466938840729e-05, /* a3 */
    1.1834955060129207e-05, /* a4 */
    5.60129697498419e-08,   /* a5 */
    -3.732041491227701e-09  /* a6 */
};
}

namespace Adsorption
{
const double DensityHauer15CaCl2CaX::M_Ads =
    MaterialLib::PhysicalConstant::MolarMass::Water;

double DensityHauer15CaCl2CaX::getAdsorbateDensity(const double T_Ads)
{
    return DensityHauer::getAdsorbateDensity(T_Ads);
}

double DensityHauer15CaCl2CaX::getAlphaT(const double T_Ads)
{
    return DensityHauer::getAlphaT(T_Ads);
}

// Characteristic curve. Return W (A)
double DensityHauer15CaCl2CaX::characteristicCurve(const double A)
{
    double W = curvePolyfrac(c, A);  // cm^3/g

    if (W < 0.0)
    {
        W = 0.0;  // TODO [CL] debug output
    }

    return W / 1.e3;  // m^3/kg
}

double DensityHauer15CaCl2CaX::dCharacteristicCurve(const double A)
{
    return dCurvePolyfrac(c, A) / 1000.0;
}
}
