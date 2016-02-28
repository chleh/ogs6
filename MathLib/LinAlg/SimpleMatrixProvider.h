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
    using MSP = MatrixSpecificationsProvider;

    Matrix& getMatrix(MSP const& msp, std::size_t& id) override;
    Matrix& getMatrix(MSP const& msp, std::size_t& id, Matrix const& A) override;

    void releaseMatrix(std::size_t const id, Matrix const& A) override;

    Vector& getVector(MSP const& msp, std::size_t& id) override;
    Vector& getVector(MSP const& msp, std::size_t& id, Vector const& x) override;
    void releaseVector(std::size_t const id, Vector const& x);

private:
    std::size_t _next_id = 1;

    std::map<std::size_t, std::unique_ptr<Matrix> > _unused_matrices;
    std::map<std::size_t, std::unique_ptr<Matrix> > _used_matrices;

    std::map<std::size_t, std::unique_ptr<Vector> > _unused_vectors;
    std::map<std::size_t, std::unique_ptr<Vector> > _used_vectors;
};



} // namespace MathLib

#include "SimpleMatrixProvider-impl.h"
