/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATERIALSLIB_ADSORPTION_DENSITYNUNEZ_H
#define MATERIALSLIB_ADSORPTION_DENSITYNUNEZ_H

#include "Adsorption.h"

namespace Adsorption
{

class DensityNunez : public AdsorptionReaction
{
public:
    DensityNunez(const double k_rate, const double p_Ads_half)
        : AdsorptionReaction(k_rate, p_Ads_half)
    {}

    double getAdsorbateDensity(const double T_Ads) const;
    double getAlphaT(const double T_Ads) const;
    double characteristicCurve(const double A) const;
    double dCharacteristicCurve(const double A) const;
};

}
#endif // MATERIALSLIB_ADSORPTION_DENSITYNUNEZ_H
