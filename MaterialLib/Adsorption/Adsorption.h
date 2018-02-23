/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <array>
#include <cmath>

namespace Adsorption
{
class AdsorptionReaction
{
public:
    static double getEquilibriumVapourPressure(const double T_Ads);
    static double getSpecificHeatCapacity(const double T_Ads); // TODO [CL] why unused?

    static double getMolarFraction(double xm, double M_this, double M_other);
    static double getMassFraction(double xn, double M_this, double M_other);
    static double dMolarFraction(double xm, double M_this, double M_other);

    static double getLoading(const double rho_curr, const double rho_dry);
};


inline double curvePolyfrac(const double* coeffs, const double x)
{
    auto const numerator =
        coeffs[0] + x * (coeffs[2] + x * (coeffs[4] + x * coeffs[6]));
    auto const denominator =
        1.0 + x * (coeffs[1] + x * (coeffs[3] + x * coeffs[5]));
    return numerator / denominator;
}

inline double dCurvePolyfrac(const double* coeffs, const double x)
{
    const double x2 = x*x;
    const double x3 = x2*x;
    const double u  = coeffs[0] + coeffs[2] * x +     coeffs[4] * x2 +     coeffs[6] * x3;
    const double du =             coeffs[2]     + 2.0*coeffs[4] * x  + 3.0*coeffs[6] * x2;
    const double v  = 1.0 + coeffs[1] * x +     coeffs[3] * x2 +     coeffs[5] * x3;
    const double dv =       coeffs[1]     + 2.0*coeffs[3] * x  + 3.0*coeffs[5] * x2;

    return (du*v - u*dv) / v / v;
}

}
