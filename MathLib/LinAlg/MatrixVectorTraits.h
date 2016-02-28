// TODO
#pragma once

#include "MatrixProviderUser.h"

namespace MathLib
{
template<typename Matrix>
struct MatrixVectorTraits;
}


// TODO add vector template param
#define SpecializeMatrixVectorTraits(MATVEC, IDX) \
    template<> struct MatrixVectorTraits<MATVEC> { \
        using Index = IDX; \
        static MATVEC* newInstance(); \
        static MATVEC* newInstance(MATVEC const& A); \
        static MATVEC* newInstance(MatrixSpecifications const& spec); \
    };


#ifdef OGS_USE_EIGEN

#include<Eigen/Core>

namespace MathLib
{
SpecializeMatrixVectorTraits(Eigen::MatrixXd, Eigen::MatrixXd::Index);
SpecializeMatrixVectorTraits(Eigen::VectorXd, Eigen::VectorXd::Index);
}

#endif


#ifdef USE_PETSC

#include "MathLib/LinAlg/PETSc/PETScMatrix.h"
#include "MathLib/LinAlg/PETSc/PETScVector.h"

namespace MathLib
{
SpecializeMatrixVectorTraits(PETScMatrix, PETScMatrix::IndexType);
SpecializeMatrixVectorTraits(PETScVector, PETScVector::IndexType);
}


#elif defined(OGS_USE_EIGEN)

#include "MathLib/LinAlg/Eigen/EigenMatrix.h"
#include "MathLib/LinAlg/Eigen/EigenVector.h"

namespace MathLib
{
SpecializeMatrixVectorTraits(EigenMatrix, EigenMatrix::IndexType);
SpecializeMatrixVectorTraits(EigenVector, EigenVector::IndexType);
}

#endif

#undef SpecializeMatrixVectorTraits
