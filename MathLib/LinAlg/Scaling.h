/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#ifdef OGS_USE_MKL
#include "MathLib/LinAlg/Sparse/LOLMatrix.h"
#include "MathLib/LinAlg/Dense/DenseVector.h"
#endif

#ifdef OGS_USE_EIGEN
#include "MathLib/LinAlg/Eigen/EigenMatrix.h"
#include "MathLib/LinAlg/Eigen/EigenVector.h"
#endif

namespace MathLib
{

#ifdef OGS_USE_MKL
void scaleDiagonal(LOLMatrix& A, DenseVector<double>& b);
#endif

#ifdef OGS_USE_EIGEN
void scaleDiagonal(EigenMatrix& A, EigenVector& b);
#endif

}
