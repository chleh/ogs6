#pragma once

#include "adsorption.h"

namespace Ads
{

class DensityDubinin : public Adsorption
{
public:
	double get_adsorbate_density(const double Tads) const;
	double get_alphaT(const double Tads) const;
	double characteristic_curve(const double A) const;
	double d_characteristic_curve(const double A) const;
};

}
