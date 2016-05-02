/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

namespace ProcessLib
{
namespace TES
{

//! Computes the dielectric heating power in W/m^3 for the given
//! temperature \c T and \c loading.
double getVolumetricJouleHeatingPower(const double T, const double loading);

} // namespace TES
} // namespace ProcessLib
