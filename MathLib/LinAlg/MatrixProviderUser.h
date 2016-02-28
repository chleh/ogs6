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

    virtual Vector& getVector() = 0;
    virtual Vector& getVector(Vector const& x) = 0;
    virtual Vector& getVector(std::size_t& id) = 0;
    //! get an uninitialized vector (or the one with the given id)
    virtual Vector& getVector(MSP const& msp, std::size_t& id) = 0;
    //! get a copy of x
    virtual Vector& getVector(MSP const& msp, std::size_t& id, Vector const& x) = 0;

    virtual void releaseVector(std::size_t const id, Vector const& x) = 0;

    virtual ~VectorProvider() = default;
};

template<typename Matrix, typename Vector>
class MatrixProvider
        : public VectorProvider<Vector>
{
public:
    using MSP = MatrixSpecificationsProvider;

    virtual Matrix& getMatrix() = 0;
    //! get an uninitialized matrix (or the one with the given id)
    virtual Matrix& getMatrix(MSP const& msp, std::size_t& id) = 0;
    //! get a copy of A
    virtual Matrix& getMatrix(MSP const& msp, std::size_t& id, Matrix const& A) = 0;

    virtual void releaseMatrix(std::size_t const id, Matrix const& A) = 0;
};

template<typename Matrix, typename Vector>
class MatrixUser
{
public:
    virtual void setMatrixProvider(MatrixProvider<Matrix, Vector>& prvd) = 0;

    virtual ~MatrixUser() = default;
};


} // namespace MathLib
