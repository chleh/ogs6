#pragma once

#include "DensityHauer.h"

namespace Adsorption
{

class DensityHauerNaYStach : public AdsorptionReaction
{
public:
    double get_adsorbate_density(const double T_Ads) const;
    double get_alphaT(const double T_Ads) const;
    double characteristic_curve(const double A) const;
    double d_characteristic_curve(const double A) const;
    double get_enthalpy(const double p_Ads, const double T_Ads, const double M_Ads) const override;
};

}
