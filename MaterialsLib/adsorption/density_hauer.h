#pragma once

#include "adsorption.h"
#include "density_cook.h"

namespace Ads
{

class DensityHauer : public Adsorption
{
public:
	double get_adsorbate_density(const double Tads) const;
	double get_alphaT(const double Tads) const;
	double characteristic_curve(const double A) const;
};

inline double rho_water_Hauer(const double Tads)
{
	// data like in python script
	const double T0 = 283.15, rho0 = rho_water_Dean(T0), alpha0 = 3.781e-4; //K; kg/m^3; 1/K

	return rho0 * (1. - alpha0 * (Tads-T0)); //in kg/m^3
}

}
