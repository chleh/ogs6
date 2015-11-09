/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LOLMatrix.h"
#include "MathLib/LinAlg/Dense/DenseVector.h"

#include <cassert>

namespace MathLib
{

void LOLMatrix::
multiply(const DenseVector<double> &x, DenseVector<double> &res) const
{
    assert(x.size() == getNCols());

    auto const nrows = getNRows();
    res.resize(nrows);

    for (IndexType r=0; r < (IndexType) nrows; ++r)
    {
        double xres = 0.0;
        auto const& row = _mat[r];

        for (auto const& e : row) {
            xres += e.value * x[e.col_idx];
        }

        res[r] = xres;
    }
}

}
