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

    Vector& getVector() override;
    Vector& getVector(std::size_t& id) override;

    Vector& getVector(Vector const& x) override;
    Vector& getVector(Vector const& x, std::size_t& id) override;

    Vector& getVector(MSP const& msp) override;
    Vector& getVector(MSP const& msp, std::size_t& id) override;

    void releaseVector(Vector const& x);

    Matrix& getMatrix() override;
    Matrix& getMatrix(std::size_t& id) override;

    Matrix& getMatrix(Matrix const& A) override;
    Matrix& getMatrix(Matrix const& A, std::size_t& id) override;

    Matrix& getMatrix(MSP const& msp) override;
    Matrix& getMatrix(MSP const& msp, std::size_t& id) override;

    void releaseMatrix(Matrix const& A) override;

    ~SimpleMatrixProvider();

private:
    template<bool do_search, typename... Args>
    Matrix& getMatrix_(std::size_t& id, Args&&... args);

    template<bool do_search, typename... Args>
    Vector& getVector_(std::size_t& id, Args&&... args);

    template<bool do_search, typename MatVec, typename... Args>
    MatVec& get_(std::size_t& id,
                 std::map<std::size_t, MatVec*>& unused_map,
                 std::map<MatVec*, std::size_t>& used_map,
                 Args&&... args);

    std::size_t _next_id = 1;

    std::map<std::size_t, Matrix*> _unused_matrices;
    std::map<Matrix*, std::size_t> _used_matrices;

    std::map<std::size_t, Vector*> _unused_vectors;
    std::map<Vector*, std::size_t> _used_vectors;
};



} // namespace MathLib

#include "SimpleMatrixProvider-impl.h"
