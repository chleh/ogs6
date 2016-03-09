/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once


#include "LOLMatrix.h"
#include "CSRMatrix.h"

namespace MathLib
{

CSRMatrix
toCSR(LOLMatrix const& mat);

LOLMatrix
toLOL(CSRMatrix const& mat);

}
