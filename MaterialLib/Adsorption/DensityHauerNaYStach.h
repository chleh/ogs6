/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATERIALLIB_ADSORPTION_DENSITYHAUERNAYSTACH_H
#define MATERIALLIB_ADSORPTION_DENSITYHAUERNAYSTACH_H

namespace Adsorption
{
class DensityHauerNaYStach
{
public:
    static double getAdsorbateDensity(const double T_Ads);
    static double getAlphaT(const double T_Ads);
    static double characteristicCurve(const double A);
    static double dCharacteristicCurve(const double A);
    static double getEnthalpy(double p_Ads, double T_Ads);
    static const double M_Ads;
};
}

#endif  // MATERIALLIB_ADSORPTION_DENSITYHAUERNAYSTACH_H
