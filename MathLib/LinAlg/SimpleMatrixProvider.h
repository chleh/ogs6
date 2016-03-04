/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATHLIB_SIMPLE_MATRIX_PROVIDER_H
#define MATHLIB_SIMPLE_MATRIX_PROVIDER_H

#include<map>
#include<memory>

#include "MatrixProviderUser.h"

namespace MathLib
{

template<typename Matrix, typename Vector>
class SimpleMatrixProvider final
        : public MatrixProvider<Matrix, Vector>
        , public VectorProvider<Vector>
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

    void releaseVector(Vector const& x) override;

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
    std::pair<Matrix*, bool> getMatrix_(std::size_t& id, Args&&... args);

    template<bool do_search, typename... Args>
    std::pair<Vector*, bool> getVector_(std::size_t& id, Args&&... args);

    // returns a pair with the pointer to the matrix/vector and
    // a boolean indicating if a new object has been built (then true else false)
    template<bool do_search, typename MatVec, typename... Args>
    std::pair<MatVec*, bool>
    get_(std::size_t& id,
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

#endif // MATHLIB_SIMPLE_MATRIX_PROVIDER_H
