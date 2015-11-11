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


LOLMatrix
LOLMatrix::transpose() const
{
    std::vector<Row> cols(_num_cols);

    std::vector<std::size_t> nnz_per_col(_num_cols);

    for (auto const& row : _mat) {
        for (auto const& e : row) {
            ++ nnz_per_col[e.col_idx];
        }
    }

    for (std::size_t c=0; c<_num_cols; ++c) {
        cols[c].reserve(nnz_per_col[c]);
    }

    auto const num_rows = _mat.size();
    for (std::size_t r=0; r<num_rows; ++r) {
        for (auto const& e : _mat[r]) {
            cols[e.col_idx].emplace_back(r, e.value);
        }
    }

    return LOLMatrix(num_rows, std::move(cols));
}

}
