/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Scaling.h"

#ifdef OGS_USE_MKL
#include <algorithm>
#endif

namespace MathLib
{

#ifdef OGS_USE_MKL
void
scaleDiagonal(LOLMatrix& A, DenseVector<double>& b)
{
    using IT = LOLMatrix::IndexType;
    for (IT r=0; r < (IT) A.getNRows(); ++r)
    {
        auto& row = A.getRow(r);
        auto d = std::lower_bound(row.cbegin(), row.cend(), r);
        if (d->col_idx == r && d->value != 0.0) {
            auto const v = d->value;
            for (auto& e : row) {
                e.value /= v;
            }
            b[r] /= v;
        }
    }
}
#endif

#ifdef OGS_USE_EIGEN
void
scaleDiagonal(EigenMatrix &A, EigenVector &b)
{
    Eigen::VectorXd diag = A.getRawMatrix().diagonal();

    for (int i=0; i<diag.size(); ++i)
    {
        A.getRawMatrix().row(i) /= diag[i];
        b.getRawVector()[i] /= diag[i];
    }
}
#endif

}
