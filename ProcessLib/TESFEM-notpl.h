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

#ifndef PROCESS_LIB_TES_FEM_NOTPL_H_
#define PROCESS_LIB_TES_FEM_NOTPL_H_

#include "TESFEM-data-notpl.h"

namespace ProcessLib
{

namespace TES
{

class LAMethodsNoTpl
{
public:
    static double getMassCoeff(LADataNoTpl const& d);
};

} // namespace TES

} // namespace ProcessLib

#endif // PROCESS_LIB_TES_FEM_NOTPL_H_
