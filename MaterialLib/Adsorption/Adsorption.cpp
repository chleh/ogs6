/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Adsorption.h"
#include <logog/include/logog.hpp>
#include "MaterialLib/PhysicalConstant.h"

namespace
{
template <typename T>
T square(const T& v)
{
    return v * v;
}
}

namespace Adsorption
{
// Saturation pressure for water used in Nunez
double AdsorptionReaction::getEquilibriumVapourPressure(const double T_Ads)
{
    // critical T and p
    const double Tc = 647.3;    // K
    const double pc = 221.2e5;  // Pa
    // dimensionless T
    const double Tr = T_Ads / Tc;
    const double theta = 1. - Tr;
    // empirical constants
    const double c[] = {-7.69123, -26.08023, -168.17065, 64.23285, -118.96462,
                        4.16717,  20.97506,  1.0e9,      6.0};
    const double K[] = {c[0] * theta + c[1] * pow(theta, 2) +
                            c[2] * pow(theta, 3) + c[3] * pow(theta, 4) +
                            c[4] * pow(theta, 5),
                        1. + c[5] * theta + c[6] * pow(theta, 2)};

    const double exponent =
        K[0] / (K[1] * Tr) - theta / (c[7] * pow(theta, 2) + c[8]);
    return pc * exp(exponent);  // in Pa
}

// evaluate specific heat capacity of adsorbate follwing Nunez
double AdsorptionReaction::getSpecificHeatCapacity(const double T_Ads)
{
    const double c[] = {4.224,         -3.716e-3,   9.351e-5,      -7.1786e-7,
                        -9.1266e-9,    2.69247e-10, -2.773104e-12, 1.553177e-14,
                        -4.982795e-17, 8.578e-20,   -6.12423e-23};
    double cp = 0.;
    for (unsigned i = 0; i < sizeof(c) / sizeof(c[0]); i++)
        cp += c[i] * pow(T_Ads, i);
    return cp;  // kJ/(kg*K)
}

double AdsorptionReaction::getMolarFraction(double xm, double M_this,
                                            double M_other)
{
    return M_other * xm / (M_other * xm + M_this * (1.0 - xm));
}

double AdsorptionReaction::getMassFraction(double xn, double M_this,
                                           double M_other)
{
    return M_this * xn / (M_this * xn + M_other * (1.0 - xn));
}

double AdsorptionReaction::dMolarFraction(double xm, double M_this,
                                          double M_other)
{
    return M_other * M_this / square(M_other * xm + M_this * (1.0 - xm));
}

double AdsorptionReaction::getLoading(const double rho_curr,
                                      const double rho_dry)
{
    return rho_curr / rho_dry - 1.0;
}

}  // namespace Ads
