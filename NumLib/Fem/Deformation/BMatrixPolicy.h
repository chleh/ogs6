/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef NUMLIB_FEM_BMATRIXPOLICY_H_
#define NUMLIB_FEM_BMATRIXPOLICY_H_

// global_dim = number of the displacement DOF
// n_nodes = nodes of an element
template <int N, int GlobalDim>
struct BMatrixDimensions
{
	static int const rows = 0;
	static int const columns = N * GlobalDim;
};

template <int N>
struct BMatrixDimensions<N, 1>
{
	static int const rows = 1;
	static int const columns = N * 1;
};

template <int N>
struct BMatrixDimensions<N, 2>
{
	static int const rows = 4;
	static int const columns = N * 2;
};

template <int N>
struct BMatrixDimensions<N, 3>
{
	static int const rows = 6;
	static int const columns = N * 3;
};

/// An implementation of ShapeMatrixPolicy using fixed size (compile-time) eigen
/// matrices and vectors.
template <typename ShapeFunction, unsigned GlobalDim>
struct EigenFixedBMatrixPolicy
{
	template <int N>
	using _VectorType = typename ::detail::EigenMatrixType<N, 1>::type;

	template <int N, int M>
	using _MatrixType = typename ::detail::EigenMatrixType<N, M>::type;

	using StiffnessMatrixType = _MatrixType<
	    BMatrixDimensions<ShapeFunction::NPOINTS, GlobalDim>::columns,
	    BMatrixDimensions<ShapeFunction::NPOINTS, GlobalDim>::columns>;
	/// rhs residual
	using NodalForceVectorType = _VectorType<
	    BMatrixDimensions<ShapeFunction::NPOINTS, GlobalDim>::columns>;

	/// sigma
	using StressVectorType =
	    _VectorType<BMatrixDimensions<ShapeFunction::NPOINTS, GlobalDim>::rows>;
	/// C
	using ModulusMatrixType =
	    _MatrixType<BMatrixDimensions<ShapeFunction::NPOINTS, GlobalDim>::rows,
	                BMatrixDimensions<ShapeFunction::NPOINTS, GlobalDim>::rows>;

	using BMatrixType = _MatrixType<
	    BMatrixDimensions<ShapeFunction::NPOINTS, GlobalDim>::rows,
	    BMatrixDimensions<ShapeFunction::NPOINTS, GlobalDim>::columns>;
};

/// An implementation of ShapeMatrixPolicy using fixed size (compile-time) eigen
/// matrices and vectors.
template <typename ShapeFunction, unsigned GlobalDim>
struct EigenDynamicBMatrixPolicy
{
	// Dynamic size local matrices are much slower in allocation than their
	// fixed counterparts.
	template <int, int>
	using _MatrixType =
	    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
	template <int>
	using _VectorType = Eigen::Matrix<double, Eigen::Dynamic, 1>;

	using StiffnessMatrixType = _MatrixType<
	    BMatrixDimensions<ShapeFunction::NPOINTS, GlobalDim>::columns,
	    BMatrixDimensions<ShapeFunction::NPOINTS, GlobalDim>::columns>;
	/// rhs residual
	using NodalForceVectorType = _VectorType<
	    BMatrixDimensions<ShapeFunction::NPOINTS, GlobalDim>::columns>;

	/// sigma
	using StressVectorType =
	    _VectorType<BMatrixDimensions<ShapeFunction::NPOINTS, GlobalDim>::rows>;
	/// C
	using ModulusMatrixType =
	    _MatrixType<BMatrixDimensions<ShapeFunction::NPOINTS, GlobalDim>::rows,
	                BMatrixDimensions<ShapeFunction::NPOINTS, GlobalDim>::rows>;

	using BMatrixType = _MatrixType<
	    BMatrixDimensions<ShapeFunction::NPOINTS, GlobalDim>::rows,
	    BMatrixDimensions<ShapeFunction::NPOINTS, GlobalDim>::columns>;
};

#ifdef OGS_EIGEN_DYNAMIC_SHAPE_MATRICES
template <typename ShapeFunction, unsigned GlobalDim>
using BMatrixPolicyType = EigenDynamicBMatrixPolicy<ShapeFunction, GlobalDim>;
#else
template <typename ShapeFunction, unsigned GlobalDim>
using BMatrixPolicyType = EigenFixedBMatrixPolicy<ShapeFunction, GlobalDim>;
#endif

#endif  // NUMLIB_FEM_BMATRIXPOLICY_H_
