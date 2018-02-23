/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "DensityHauerNaYStachCubicAndExp.h"
#include <cmath>
#include "DensityHauer.h"
#include "MaterialLib/PhysicalConstant.h"

// out-TCP/NaY_Stach_Hauer_cubicAndExp_CC.pickle
// date extracted 2016-07-12 15:47:39 file mtime 2016-07-12 15:01:27
const double CC_A0 = 670.4032522615867;
const double CC_a3 = -1.2959841712540927e-10;
const double CC_b = -445.0325080309215;
const double CC_bp = 0.0029233320535917776;

const double CC_A0_m_b = CC_A0 - CC_b;
const double CC_A0_m_b2 = CC_A0_m_b * CC_A0_m_b;
const double CC_a0 = (-3.0 / CC_bp - CC_A0_m_b) * CC_a3 * CC_A0_m_b2;
const double CC_ap =
    3.0 * CC_a0 * std::exp(CC_bp * CC_A0) / (3.0 + CC_bp * CC_A0_m_b);

template <int i>
double mypow(const double x)
{
    if (i < 0)
    {
        return 1.0 / mypow<-i>(x);
    }
    else
    {
        const double p = mypow<(i >> 1)>(x);
        return (i & 1) ? p * p * x : p * p;
    }
}

template <>
inline double mypow<0>(const double /*x*/)
{
    return 1.0;
}

namespace Adsorption
{
const double DensityHauerNaYStachCubicAndExp::M_Ads =
    MaterialLib::PhysicalConstant::MolarMass::Water;

double DensityHauerNaYStachCubicAndExp::getAdsorbateDensity(const double T_Ads)
{
    return DensityHauer::getAdsorbateDensity(T_Ads);
}

// Thermal expansivity model for water found in the works of Hauer
double DensityHauerNaYStachCubicAndExp::getAlphaT(const double T_Ads)
{
    return DensityHauer::getAlphaT(T_Ads);
}

// Characteristic curve. Return W (A)
double DensityHauerNaYStachCubicAndExp::characteristicCurve(const double A)
{
    double W;
    if (A <= CC_A0)
    {
        W = CC_a3 * mypow<3>(A - CC_b) + CC_a0;
    }
    else
    {
        W = CC_ap * std::exp(-CC_bp * A);
    }

    return W / 1.e3;  // m^3/kg
}

double DensityHauerNaYStachCubicAndExp::dCharacteristicCurve(const double A)
{
    double W;
    if (A <= CC_A0)
    {
        W = 3 * CC_a3 * mypow<2>(A - CC_b);
    }
    else
    {
        W = -CC_ap * CC_bp * std::exp(-CC_bp * A);
    }
    return W / 1000.0;
}

double DensityHauerNaYStachCubicAndExp::getEnthalpy(double, double)
{
    return 2055.0 * 1000.0;  // J/kg value taken from Kraus' PhD thesis
}
}
