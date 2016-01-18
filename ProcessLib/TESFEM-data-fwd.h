/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_TESFEM_DATA_FWD_H_
#define PROCESS_LIB_TESFEM_DATA_FWD_H_

#include "TESFEM-data.h"

namespace ProcessLib
{

namespace TES
{

#ifdef EIGEN_DYNAMIC_SHAPE_MATRICES
extern template class LADataNoTpl<DataTraits<int, 0, 0, 0> >;
static_assert(EIGEN_DYNAMIC_SHAPE_MATRICES_FLAG == 1, "inconsistent use of macros");
#else
static_assert(EIGEN_DYNAMIC_SHAPE_MATRICES_FLAG == 0, "inconsistent use of macros");
#endif

}

}

#endif  // PROCESS_LIB_TESFEM_DATA_FWD_H_
