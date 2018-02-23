#pragma once

#include "DensityHauer.h"

namespace Adsorption
{
class DensityHauerCaX80NoOutliers
{
public:
    static double getAdsorbateDensity(const double T_Ads);
    static double getAlphaT(const double T_Ads);
    static double characteristicCurve(const double A);
    static double dCharacteristicCurve(const double A);
    static const double M_Ads;
};
}
