// TODO
#pragma once

#include <cstddef>

namespace MathLib
{

struct MatrixSpecifications
{
    std::size_t const nrows;
    std::size_t const ncols;
};

class MatrixSpecificationsProvider
{
public:
    virtual MatrixSpecifications getMatrixSpecifications() const = 0;

    virtual ~MatrixSpecificationsProvider() = default;
};


template<typename Vector>
class VectorProvider
{
public:
    using MSP = MatrixSpecificationsProvider;

    //! get an uninitialized vector (or the one with the given id)
    virtual Vector& getVector() = 0;
    virtual Vector& getVector(std::size_t& id) = 0;

    //! get a copy of x
    virtual Vector& getVector(Vector const& x) = 0;
    virtual Vector& getVector(Vector const& x, std::size_t& id) = 0;

    //! get a vector according to the given specifications
    virtual Vector& getVector(MSP const& msp) = 0;
    virtual Vector& getVector(MSP const& msp, std::size_t& id) = 0;

    virtual void releaseVector(Vector const& x) = 0;

    virtual ~VectorProvider() = default;
};

template<typename Matrix, typename Vector>
class MatrixProvider
{
public:
    using MSP = MatrixSpecificationsProvider;

    //! get an uninitialized Matrix (or the one with the given id)
    virtual Matrix& getMatrix() = 0;
    virtual Matrix& getMatrix(std::size_t& id) = 0;

    //! get a copy of x
    virtual Matrix& getMatrix(Matrix const& A) = 0;
    virtual Matrix& getMatrix(Matrix const& A, std::size_t& id) = 0;

    //! get a Matrix according to the given specifications
    virtual Matrix& getMatrix(MSP const& msp) = 0;
    virtual Matrix& getMatrix(MSP const& msp, std::size_t& id) = 0;

    virtual void releaseMatrix(Matrix const& A) = 0;
};

} // namespace MathLib
