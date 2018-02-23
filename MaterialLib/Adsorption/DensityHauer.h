/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "DensityCook.h"

namespace Adsorption
{
struct DensityHauer
{
    static double getAdsorbateDensity(const double T_Ads);
    static double getAlphaT(const double T_Ads);
    static double characteristicCurve(const double A);
    static double dCharacteristicCurve(const double A);
    static const double M_Ads;
};
}
