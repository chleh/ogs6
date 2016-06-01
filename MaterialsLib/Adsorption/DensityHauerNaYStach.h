#pragma once

#include "DensityHauer.h"

namespace Adsorption
{

class DensityHauerNaYStach : public AdsorptionReaction
{
public:
    double getAdsorbateDensity(const double T_Ads) const;
    double getAlphaT(const double T_Ads) const;
    double characteristicCurve(const double A) const;
    double dCharacteristicCurve(const double A) const;
    double getEnthalpy(const double p_Ads, const double T_Ads, const double M_Ads) const override;
};

}
