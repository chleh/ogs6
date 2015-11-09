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


namespace MathLib
{

void applyKnownSolution(LOLMatrix &A, DenseVector<double> &b,
                        const std::vector<std::size_t> &vec_knownX_id,
                        const std::vector<double> &vec_knownX_x,
                        double /*penalty_scaling*/)
{
    using IT = LOLMatrix::IndexType;

    auto const n_rows = A.getNRows();

    for (std::size_t ix=0; ix<vec_knownX_id.size(); ix++)
    {
        IT const row_id = vec_knownX_id[ix];
        auto const x = vec_knownX_x[ix];

        //A(k, j) = 0.
        for (auto& e : A.getRow(row_id)) e.value = 0.0;
        // A.setRowZero(row_id);

        //b_i -= A(i,k)*val, i!=k
        for (std::size_t i=0; i<n_rows; i++)
        {
            auto& row = A.getRow(i);

            for (auto& e : row) {
                if (e.col_idx == row_id) {
                    b[i] -= e.value * x;
                    e.value = 0.0;
                    break;
                }
            }
        }

        //A(k, k) = 1.0
        A.setValue(row_id, row_id, 1.0);

        //b_k = val
        b[row_id] = x;
    }
}

}
