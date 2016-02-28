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

template<typename Matrix, typename Vector>
class MatrixProvider
{
public:
    virtual Matrix& getMatrix(std::size_t const size, std::size_t& id) = 0;
    virtual void releaseMatrix(std::size_t const id) = 0;

    virtual Vector& getVector(std::size_t const size, std::size_t& id) = 0;
    virtual void releaseVector(std::size_t const id) = 0;

    virtual void setMatrixSpecificationsProvider(MatrixSpecificationsProvider const& spec_prvd) = 0;

    virtual ~MatrixProvider() = default;
};

template<typename Matrix, typename Vector>
class MatrixUser
{
public:
    virtual void setMatrixProvider(MatrixProvider<Matrix, Vector>& prvd) = 0;

    virtual ~MatrixUser() = default;
};


} // namespace MathLib
