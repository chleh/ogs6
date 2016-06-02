/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BLAS.h"

namespace MathLib
{
namespace BLAS
{
namespace detail
{
template <typename MatVec>
class BLASHelper;
}  // namespace detail
}  // namespace BLAS
}  // namespace MathLib


// TODO reorder BLAS function signatures?

// Dense Eigen matrices and vectors ////////////////////////////////////////
#ifdef OGS_USE_EIGEN

#include <Eigen/Core>

namespace MathLib
{
namespace BLAS
{

template <>
void setZero(Eigen::VectorXd& x)
{
    x.setZero();
}

template <>
void setZero(Eigen::MatrixXd& A)
{
    A.setZero();
}

// Explicit specialization
// Computes w = x/y componentwise.
template<>
void componentwiseDivide(Eigen::VectorXd& w,
                         Eigen::VectorXd const& x, Eigen::VectorXd const& y)
{
    w.noalias() = x.cwiseQuotient(y);
}

// Explicit specialization
// Computes the Manhattan norm of x
template<>
double norm1(Eigen::VectorXd const& x)
{
    return x.lpNorm<1>();
}

// Explicit specialization
// Computes the Euclidean norm of x
template<>
double norm2(Eigen::VectorXd const& x)
{
    return x.norm();
}

// Explicit specialization
// Computes the Maximum norm of x
template<>
double normMax(Eigen::VectorXd const& x)
{
    return x.lpNorm<Eigen::Infinity>();
}

template<>
MatrixVectorTraits<Eigen::VectorXd>::Index sizeGlobal(Eigen::VectorXd const& x)
{
    return x.size();
}
} } // namespaces

#endif


// Global PETScMatrix/PETScVector //////////////////////////////////////////
#ifdef USE_PETSC

#include "MathLib/LinAlg/PETSc/PETScVector.h"
#include "MathLib/LinAlg/PETSc/PETScMatrix.h"

namespace MathLib
{
namespace BLAS
{
namespace detail
{
template <>
class BLASHelper<PETScVector>
{
public:
    static MatrixVectorTraits<PETScVector>::Index sizeGlobal(
        PETScVector const& x)
    {
        return x._size;
    }

    static MatrixVectorTraits<PETScVector>::Index sizeLocalWithoutGhosts(
        PETScVector const& x)
    {
        return x._size_loc;
    }

    static MatrixVectorTraits<PETScVector>::Index numberOfGhosts(
        PETScVector const& x)
    {
        return x._size_ghosts;
    }

    static void copy(PETScVector const& x, PETScVector& y)
    {
        if (!y._v) y.shallowCopy(x);
        VecCopy(*x._v, *y._v);
    }

    static void setZero(PETScVector& x)
    {
        VecSet(*x._v, 0.0);
    }
};
}  // namespace detail

// Vector

template <>
void setZero(PETScVector& x)
{
    detail::BLASHelper<PETScVector>::setZero(x);
}

template <>
void setZero(PETScMatrix& A)
{
    MatZeroEntries(A.getRawMatrix());
}

void copy(PETScVector const& x, PETScVector& y)
{
    detail::BLASHelper<PETScVector>::copy(x, y);
}

void scale(PETScVector& x, double const a)
{
    VecScale(*x.getRawVector(), a);
}

// y = a*y + X
void aypx(PETScVector& y, double const a, PETScVector const& x)
{
    // TODO check sizes
    VecAYPX(*y.getRawVector(), a, *x.getRawVector());
}

// y = a*x + y
void axpy(PETScVector& y, double const a, PETScVector const& x)
{
    // TODO check sizes
    VecAXPY(*y.getRawVector(), a, *x.getRawVector());
}

// y = a*x + b*y
void axpby(PETScVector& y, double const a, double const b, PETScVector const& x)
{
    // TODO check sizes
    VecAXPBY(*y.getRawVector(), a, b, *x.getRawVector());
}

// Explicit specialization
// Computes w = x/y componentwise.
template<>
void componentwiseDivide(PETScVector& w,
                         PETScVector const& x, PETScVector const& y)
{
    VecPointwiseDivide(*w.getRawVector(), *x.getRawVector(), *y.getRawVector());
}

// Explicit specialization
// Computes the Manhattan norm of x
template<>
double norm1(PETScVector const& x)
{
    PetscScalar norm = 0;
    VecNorm(*x.getRawVector(), NORM_1, &norm);
    return norm;
}

// Explicit specialization
// Computes the Euclidean norm of x
template<>
double norm2(PETScVector const& x)
{
    PetscScalar norm = 0;
    VecNorm(*x.getRawVector(), NORM_2, &norm);
    return norm;
}

// Explicit specialization
// Computes the Maximum norm of x
template<>
double normMax(PETScVector const& x)
{
    PetscScalar norm = 0;
    VecNorm(*x.getRawVector(), NORM_INFINITY, &norm);
    return norm;
}


// Matrix

void copy(PETScMatrix const& A, PETScMatrix& B)
{
    B = A;
}

// A = a*A
void scale(PETScMatrix& A, double const a)
{
    MatScale(A.getRawMatrix(), a);
}

// Y = a*Y + X
void aypx(PETScMatrix& Y, double const a, PETScMatrix const& X)
{
    // TODO check sizes
    // TODO sparsity pattern, currently they are assumed to be different (slow)
    MatAYPX(Y.getRawMatrix(), a, X.getRawMatrix(),
            DIFFERENT_NONZERO_PATTERN);
}

// Y = a*X + Y
void axpy(PETScMatrix& Y, double const a, PETScMatrix const& X)
{
    // TODO check sizes
    // TODO sparsity pattern, currently they are assumed to be different (slow)
    MatAXPY(Y.getRawMatrix(), a, X.getRawMatrix(),
            DIFFERENT_NONZERO_PATTERN);
}


// Matrix and Vector

// v3 = A*v1 + v2
void matMult(PETScMatrix const& A, PETScVector const& x, PETScVector& y)
{
    // TODO check sizes
    assert(&x != &y);
    if (!y.getRawVector()) y.shallowCopy(x);
    MatMult(A.getRawMatrix(), *x.getRawVector(), *y.getRawVector());
}

// v3 = A*v1 + v2
void matMultAdd(PETScMatrix const& A, PETScVector const& v1,
                       PETScVector const& v2, PETScVector& v3)
{
    // TODO check sizes
    assert(&v1 != &v3);
    if (!v3.getRawVector()) v3.shallowCopy(v1);
    MatMultAdd(A.getRawMatrix(), *v1.getRawVector(), *v2.getRawVector(), *v3.getRawVector());
}

void finalizeAssembly(PETScMatrix& A)
{
    A.finalizeAssembly(MAT_FINAL_ASSEMBLY);
}

void finalizeAssembly(PETScVector& x)
{
    x.finalizeAssembly();
}

template<>
MatrixVectorTraits<PETScVector>::Index sizeGlobal(PETScVector const& x)
{
    return detail::BLASHelper<PETScVector>::sizeGlobal(x);
}

MatrixVectorTraits<PETScVector>::Index sizeLocalWithoutGhosts(
    PETScVector const& x)
{
    return detail::BLASHelper<PETScVector>::sizeLocalWithoutGhosts(x);
}

MatrixVectorTraits<PETScVector>::Index sizeLocalWithGhosts(PETScVector const& x)
{
    return detail::BLASHelper<PETScVector>::sizeLocalWithoutGhosts(x)
            + detail::BLASHelper<PETScVector>::numberOfGhosts(x);
}

MatrixVectorTraits<PETScVector>::Index numberOfGhosts(PETScVector const& x)
{
    return detail::BLASHelper<PETScVector>::numberOfGhosts(x);;
}

double getComponent(PETScVector const& x,
                    MatrixVectorTraits<PETScVector>::Index const i)
{
    PetscScalar value;
    VecGetValues(*x.getRawVector(), 1, &i, &value);
    return value;
}
}} // namespaces



// Sparse global EigenMatrix/EigenVector //////////////////////////////////////////
#elif defined(OGS_USE_EIGEN)

#include "MathLib/LinAlg/Eigen/EigenVector.h"
#include "MathLib/LinAlg/Eigen/EigenMatrix.h"

namespace MathLib { namespace BLAS
{

template <>
void setZero(EigenVector& x)
{
    x.getRawVector().setZero();
}

template <>
void setZero(EigenMatrix& A)
{
    A.setZero();
}

// Vector

void copy(EigenVector const& x, EigenVector& y)
{
    y.getRawVector() = x.getRawVector();
}

void scale(EigenVector& x, double const a)
{
    x.getRawVector() *= a;
}

// y = a*y + X
void aypx(EigenVector& y, double const a, EigenVector const& x)
{
    // TODO: does that break anything?
    y.getRawVector() = a*y.getRawVector() + x.getRawVector();
}

// y = a*x + y
void axpy(EigenVector& y, double const a, EigenVector const& x)
{
    // TODO: does that break anything?
    y.getRawVector() += a*x.getRawVector();
}

// y = a*x + y
void axpby(EigenVector& y, double const a, double const b, EigenVector const& x)
{
    // TODO: does that break anything?
    y.getRawVector() = a*x.getRawVector() + b*y.getRawVector();
}

// Explicit specialization
// Computes w = x/y componentwise.
template<>
void componentwiseDivide(EigenVector& w,
                         EigenVector const& x, EigenVector const& y)
{
    w.getRawVector().noalias() =
            x.getRawVector().cwiseQuotient(y.getRawVector());
}

// Explicit specialization
// Computes the Manhattan norm of x
template<>
double norm1(EigenVector const& x)
{
    return x.getRawVector().lpNorm<1>();
}

// Explicit specialization
// Euclidean norm
template<>
double norm2(EigenVector const& x)
{
    return x.getRawVector().norm();
}

// Explicit specialization
// Computes the Maximum norm of x
template<>
double normMax(EigenVector const& x)
{
    return x.getRawVector().lpNorm<Eigen::Infinity>();
}


// Matrix

void copy(EigenMatrix const& A, EigenMatrix& B)
{
    B = A;
}

// A = a*A
void scale(EigenMatrix& A, double const a)
{
    // TODO: does that break anything?
    A.getRawMatrix() *= a;
}

// Y = a*Y + X
void aypx(EigenMatrix& Y, double const a, EigenMatrix const& X)
{
    // TODO: does that break anything?
    Y.getRawMatrix() = a*Y.getRawMatrix() + X.getRawMatrix();
}

// Y = a*X + Y
void axpy(EigenMatrix& Y, double const a, EigenMatrix const& X)
{
    // TODO: does that break anything?
    Y.getRawMatrix() = a*X.getRawMatrix() + Y.getRawMatrix();
}


// Matrix and Vector

// v3 = A*v1 + v2
void matMult(EigenMatrix const& A, EigenVector const& x, EigenVector& y)
{
    assert(&x != &y);
    A.multiply(x, y);
}

// v3 = A*v1 + v2
void matMultAdd(EigenMatrix const& A, EigenVector const& v1, EigenVector const& v2, EigenVector& v3)
{
    assert(&v1 != &v3);
    // TODO: does that break anything?
    v3.getRawVector() = v2.getRawVector() + A.getRawMatrix()*v1.getRawVector();
}

void finalizeAssembly(EigenMatrix& x)
{
    x.getRawMatrix().makeCompressed();
}

double getComponent(EigenVector const& x,
                    MatrixVectorTraits<EigenVector>::Index const i)
{
    return x.getRawVector()[i];
}

template<>
MatrixVectorTraits<EigenVector>::Index sizeGlobal(EigenVector const& x)
{
    return x.getRawVector().size();
}

} // namespace BLAS

} // namespace MathLib

#endif
