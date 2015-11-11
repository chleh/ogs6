/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "PardisoTools.h"
#include "MathLib/LinAlg/Sparse/LOLMatrix.h"

#include "logog/include/logog.hpp"

namespace MathLib
{

void applyKnownSolution(LOLMatrix &A, DenseVector<double> &b,
                        const std::vector<std::size_t> &vec_knownX_id,
                        const std::vector<double> &vec_knownX_x,
                        double /*penalty_scaling*/)
{
    using IT = LOLMatrix::IndexType;

    //A(k, j) = 0.
    // set row to zero
    // A.setRowZero(row_id);
    for (auto row_id : vec_knownX_id)
        for (auto& e : A.getRow(row_id)) e.value = 0.0;

    auto AT = A.transpose();

    for (std::size_t ix=0; ix<vec_knownX_id.size(); ix++)
    {
        IT const row_id = vec_knownX_id[ix];
        auto const x = vec_knownX_x[ix];

        // b_i -= A(i,k)*val, i!=k
        // set column to zero, subtract from rhs
        for (auto& e : AT.getRow(row_id)) {
            b[e.col_idx] -= e.value * x;
            e.value = 0.0;
        }

        b[row_id] = x;
        AT.setValue(row_id, row_id, 1.0);
    }

    A = AT.transpose();
}

}
