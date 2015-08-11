/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * The code of this file is used to decouple the evaluation of matrix elements from the rest of OGS6,
 * not all of OGS6 has to be recompiled every time a small change is done.
 */

#include "TESFEM-notpl.h"

namespace ProcessLib
{

namespace TES
{

double
LAMethodsNoTpl::
getMassCoeff(LADataNoTpl const& d)
{
    return d._hydraulic_conductivity;
}

} // namespace TES

} // namespace ProcessLib

