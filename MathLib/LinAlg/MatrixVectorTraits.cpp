// TODO docs


#include "MatrixVectorTraits.h"

#ifdef OGS_USE_EIGEN

namespace MathLib
{

Eigen::MatrixXd*
MatrixVectorTraits<Eigen::MatrixXd>::
newInstance()
{
    return new Eigen::MatrixXd;
}

Eigen::MatrixXd*
MatrixVectorTraits<Eigen::MatrixXd>::
newInstance(Eigen::MatrixXd const& A)
{
    return new Eigen::MatrixXd(A);
}

Eigen::MatrixXd*
MatrixVectorTraits<Eigen::MatrixXd>::
newInstance(MatrixSpecifications const& spec)
{
    return new Eigen::MatrixXd(spec.nrows, spec.ncols);
}

Eigen::VectorXd*
MatrixVectorTraits<Eigen::VectorXd>::
newInstance()
{
    return new Eigen::VectorXd;
}

Eigen::VectorXd*
MatrixVectorTraits<Eigen::VectorXd>::
newInstance(Eigen::VectorXd const& A)
{
    return new Eigen::VectorXd(A);
}

Eigen::VectorXd*
MatrixVectorTraits<Eigen::VectorXd>::
newInstance(MatrixSpecifications const& spec)
{
    return new Eigen::VectorXd(spec.nrows);
}

} // namespace MathLib

#endif // OGS_USE_EIGEN


#ifdef USE_PETSC

namespace MathLib
{

PETScMatrix*
MatrixVectorTraits<PETScMatrix>::
newInstance()
{
    return new PETScMatrix(0, 0); // TODO default constructor
}

PETScMatrix*
MatrixVectorTraits<PETScMatrix>::
newInstance(PETScMatrix const& A)
{
    return new PETScMatrix(A);
}

PETScMatrix*
MatrixVectorTraits<PETScMatrix>::
newInstance(MatrixSpecifications const& spec)
{
    return new PETScMatrix(spec.nrows); // TODO sparsity pattern
}

PETScVector*
MatrixVectorTraits<PETScVector>::
newInstance()
{
    return new PETScVector;
}

PETScVector*
MatrixVectorTraits<PETScVector>::
newInstance(PETScVector const& x)
{
    return new PETScVector(x);
}

PETScVector*
MatrixVectorTraits<PETScVector>::
newInstance(MatrixSpecifications const& spec)
{
    return new PETScVector(spec.nrows);
}

} // namespace MathLib


#elif defined(OGS_USE_EIGEN)

// #include "MathLib/LinAlg/Eigen/EigenMatrix.h"
// #include "MathLib/LinAlg/Eigen/EigenVector.h"

namespace MathLib
{

EigenMatrix*
MatrixVectorTraits<EigenMatrix>::
newInstance()
{
    return new EigenMatrix(0, 0); // TODO default constructor
}

EigenMatrix*
MatrixVectorTraits<EigenMatrix>::
newInstance(EigenMatrix const& A)
{
    return new EigenMatrix(A);
}

EigenMatrix*
MatrixVectorTraits<EigenMatrix>::
newInstance(MatrixSpecifications const& spec)
{
    return new EigenMatrix(spec.nrows); // TODO sparsity pattern
}

EigenVector*
MatrixVectorTraits<EigenVector>::
newInstance()
{
    return new EigenVector;
}

EigenVector*
MatrixVectorTraits<EigenVector>::
newInstance(EigenVector const& x)
{
    return new EigenVector(x);
}

EigenVector*
MatrixVectorTraits<EigenVector>::
newInstance(MatrixSpecifications const& spec)
{
    return new EigenVector(spec.nrows);
}

} // namespace MathLib

#endif // defined(OGS_USE_EIGEN)
