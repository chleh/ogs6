#pragma once

namespace Adsorption
{
class DensityHauer15CaCl2CaX
{
public:
    static double getAdsorbateDensity(const double T_Ads);
    static double getAlphaT(const double T_Ads);
    static double characteristicCurve(const double A);
    static double dCharacteristicCurve(const double A);
    static const double M_Ads;
};
}
