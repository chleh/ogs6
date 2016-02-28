// TODO
#pragma once

namespace MathLib
{
template<typename Matrix>
struct MatrixTraits
/*
{
    using Index = int;
}
// */
;
}


#ifdef OGS_USE_EIGEN

#include<Eigen/Core>

namespace MathLib
{
template<>
struct MatrixTraits<Eigen::MatrixXd>
{
    using Index = Eigen::MatrixXd::Index;
};
}

#endif


#ifdef USE_PETSC

namespace MathLib
{
class PETScMatrix;

template<>
struct MatrixTraits<MathLib::PETScMatrix>
{
    using Index = MathLib::PETScMatrix::IndexType;
};
}

#elif defined(OGS_USE_EIGEN)

namespace MathLib
{
class EigenMatrix;

template<>
struct MatrixTraits<MathLib::EigenMatrix>
{
    using Index = MathLib::EigenMatrix::IndexType;
};
}

#endif
