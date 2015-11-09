/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "PardisoTools.h"


namespace MathLib
{

void applyKnownSolution(LOLMatrix &A, DenseVector<double> &b,
                        const std::vector<std::size_t> &_vec_knownX_id,
                        const std::vector<double> &_vec_knownX_x,
                        double /*penalty_scaling*/)
{
    // TODO implement

}

}
