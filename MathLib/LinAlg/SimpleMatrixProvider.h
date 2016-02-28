// TODO
#pragma once

#include<map>
#include<memory>

#include "MatrixProviderUser.h"

namespace MathLib
{

template<typename Matrix, typename Vector>
class SimpleMatrixProvider final
        : public MatrixProvider<Matrix, Vector>
{
public:
    SimpleMatrixProvider() = default;

    // no copies
    SimpleMatrixProvider(SimpleMatrixProvider const&) = delete;
    SimpleMatrixProvider& operator=(SimpleMatrixProvider const&) = delete;

    using MSP = MatrixSpecificationsProvider;

    Matrix& getMatrix() override;
    Matrix& getMatrix(MSP const& msp, std::size_t& id) override;
    Matrix& getMatrix(MSP const& msp, std::size_t& id, Matrix const& A) override;

    void releaseMatrix(std::size_t const id, Matrix const& A) override;

    Vector& getVector() override;
    Vector& getVector(Vector const& x) override;
    Vector& getVector(std::size_t& id) override;
    Vector& getVector(MSP const& msp, std::size_t& id) override;
    Vector& getVector(MSP const& msp, std::size_t& id, Vector const& x) override;
    void releaseVector(std::size_t const id, Vector const& x);

    ~SimpleMatrixProvider();

private:
    std::size_t _next_id = 1;

    std::map<std::size_t, Matrix*> _unused_matrices;
    std::map<Matrix*, std::size_t> _used_matrices;

    std::map<std::size_t, Vector*> _unused_vectors;
    std::map<Vector*, std::size_t> _used_vectors;
};



} // namespace MathLib

#include "SimpleMatrixProvider-impl.h"
