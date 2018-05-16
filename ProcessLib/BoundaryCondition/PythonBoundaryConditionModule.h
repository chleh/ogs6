/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#if OGS_USE_PYTHON

#include <pybind11/pybind11.h>

namespace ProcessLib
{
void pythonBindBoundaryCondition(pybind11::module& m);

}  // namespace ProcessLib

#endif
