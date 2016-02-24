/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef NUMLIB_FEM_LINEARBMATRIX_H_
#define NUMLIB_FEM_LINEARBMATRIX_H_

namespace LinearBMatrix
{

/// Fills a B-matrix based on given shape function dN/dx values.
template <int DisplacementDim, int NPOINTS, typename GlobalDimNodalMatrixType,
          typename BMatrixType>
void computeBMatrix(GlobalDimNodalMatrixType const& dNdx, BMatrixType& b_matrix)
{
	assert(0 < DisplacementDim && DisplacementDim <= 3);

	b_matrix.setZero();

	switch (DisplacementDim)
	{
		case 3:
		{
			for (int i = 0; i < NPOINTS; ++i)
			{
				b_matrix(2, 2 * NPOINTS + i) = dNdx(2, i);
				b_matrix(4,     NPOINTS + i) = dNdx(2, i) / std::sqrt(2);
				b_matrix(4, 2 * NPOINTS + i) = dNdx(1, i) / std::sqrt(2);
				b_matrix(5,               i) = dNdx(2, i) / std::sqrt(2);
				b_matrix(5, 2 * NPOINTS + i) = dNdx(0, i) / std::sqrt(2);
			}
		}
		case 2:
		{
			for (int i = 0; i < NPOINTS; ++i)
			{
				b_matrix(1, NPOINTS + i) = dNdx(1, i);
				b_matrix(3,           i) = dNdx(1, i) / std::sqrt(2);
				b_matrix(3, NPOINTS + i) = dNdx(0, i) / std::sqrt(2);
			}
		}
		// no break for fallthrough.
		case 1:
		{
			for (int i = 0; i < NPOINTS; ++i)
				b_matrix(0, i) = dNdx(0, i);
			break;
		}
		default:
			break;
	}
}

template <typename NodalDisplacementVectorType, typename BMatrixType,
          typename StrainVectorType>
void computeStrain(NodalDisplacementVectorType const& local_displacement,
                         BMatrixType const& b_matrix, StrainVectorType& eps)
{
	eps.noalias() = b_matrix * local_displacement;
}

}   // namespace
#endif  // NUMLIB_FEM_LINEARBMATRIX_H_
