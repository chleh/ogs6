// TODO doc

#include <cassert>
#include <logog/include/logog.hpp>

#include "BLAS.h"
#include "MatrixVectorTraits.h"
#include "SimpleMatrixProvider.h"

namespace detail
{

template<typename MatVec>
MatVec*
transfer(std::map<std::size_t, MatVec*>& from_unused,
         std::map<MatVec*, std::size_t>& to_used,
         typename std::map<std::size_t, MatVec*>::iterator it)
{
    auto const id = it->first;
    auto& ptr = it->second;

    auto res = to_used.emplace(std::move(ptr), id);
    assert(res.second && "Emplacement failed.");
    from_unused.erase(it);
    return res.first->first;
}

template<typename MatVec>
void
transfer(std::map<MatVec*, std::size_t>& from_used,
         std::map<std::size_t, MatVec*>& to_unused,
         typename std::map<MatVec*, std::size_t>::iterator it)
{
    auto& ptr = it->first;
    auto const id = it->second;

    auto res = to_unused.emplace(id, std::move(ptr));
    assert(res.second && "Emplacement failed.");
    from_used.erase(it);
}

} // detail


namespace MathLib
{

template<typename Matrix, typename Vector>
template<bool do_search, typename MatVec, typename... Args>
std::pair<MatVec*, bool>
SimpleMatrixProvider<Matrix, Vector>::
get_(std::size_t& id,
     std::map<std::size_t, MatVec*>& unused_map,
     std::map<MatVec*, std::size_t>& used_map,
     Args&&... args)
{
    if (do_search)
    {
        auto it = unused_map.find(id);
        if (it != unused_map.end()) // unused matrix/vector found
            return { ::detail::transfer(unused_map, used_map, it), false };
    }

    // not searched or not found, so create a new one
    id = _next_id++;
    auto res = used_map.emplace(
        MatrixVectorTraits<MatVec>::newInstance(std::forward<Args>(args)...), id);
    assert(res.second && "Emplacement failed.");
    return { res.first->first, true };
}

template<typename Matrix, typename Vector>
template<bool do_search, typename... Args>
std::pair<Matrix*, bool>
SimpleMatrixProvider<Matrix, Vector>::
getMatrix_(std::size_t& id, Args&&... args)
{
    return get_<do_search>(id, _unused_matrices, _used_matrices, std::forward<Args>(args)...);
}


template<typename Matrix, typename Vector>
Matrix&
SimpleMatrixProvider<Matrix, Vector>::
getMatrix()
{
    std::size_t id;
    return *getMatrix_<false>(id).first;
}

template<typename Matrix, typename Vector>
Matrix&
SimpleMatrixProvider<Matrix, Vector>::
getMatrix(std::size_t& id)
{
    return *getMatrix_<true>(id).first;
}

template<typename Matrix, typename Vector>
Matrix&
SimpleMatrixProvider<Matrix, Vector>::
getMatrix(MatrixSpecificationsProvider const& msp)
{
    std::size_t id;
    auto const mat_spec = msp.getMatrixSpecifications();
    return *getMatrix_<false>(id, mat_spec).first;
    // TODO assert that the returned object always is of the right size
}

template<typename Matrix, typename Vector>
Matrix&
SimpleMatrixProvider<Matrix, Vector>::
getMatrix(MatrixSpecificationsProvider const& msp, std::size_t& id)
{
    auto mat_spec = msp.getMatrixSpecifications();
    return *getMatrix_<true>(id, mat_spec).first;
    // TODO assert that the returned object always is of the right size
}

template<typename Matrix, typename Vector>
Matrix&
SimpleMatrixProvider<Matrix, Vector>::
getMatrix(Matrix const& A)
{
    std::size_t id;
    auto const& res = getMatrix_<false>(id, A);
    if (!res.second) // no new object has been created
        BLAS::copy(A, *res.first);
    return *res.first;
}

template<typename Matrix, typename Vector>
Matrix&
SimpleMatrixProvider<Matrix, Vector>::
getMatrix(Matrix const& A, std::size_t& id)
{
    auto const& res = getMatrix_<false>(id, A);
    if (!res.second) // no new object has been created
        BLAS::copy(A, *res.first);
    return *res.first;
}

template<typename Matrix, typename Vector>
void
SimpleMatrixProvider<Matrix, Vector>::
releaseMatrix(Matrix const& A)
{
    auto it = _used_matrices.find(const_cast<Matrix*>(&A));
    if (it == _used_matrices.end()) {
        ERR("A matrix with the id %lu has not been found. Cannot release it. Aborting.");
        std::abort();
    } else {
        ::detail::transfer(_used_matrices, _unused_matrices, it);
    }
}

template<typename Matrix, typename Vector>
template<bool do_search, typename... Args>
std::pair<Vector*, bool>
SimpleMatrixProvider<Matrix, Vector>::
getVector_(std::size_t& id, Args&&... args)
{
    return get_<do_search>(id, _unused_vectors, _used_vectors, std::forward<Args>(args)...);
}


template<typename Matrix, typename Vector>
Vector&
SimpleMatrixProvider<Matrix, Vector>::
getVector()
{
    std::size_t id;
    return *getVector_<false>(id).first;
}

template<typename Matrix, typename Vector>
Vector&
SimpleMatrixProvider<Matrix, Vector>::
getVector(std::size_t& id)
{
    return *getVector_<true>(id).first;
}

template<typename Matrix, typename Vector>
Vector&
SimpleMatrixProvider<Matrix, Vector>::
getVector(MatrixSpecificationsProvider const& msp)
{
    std::size_t id;
    auto const mat_spec = msp.getMatrixSpecifications();
    return *getVector_<false>(id, mat_spec).first;
    // TODO assert that the returned object always is of the right size
}

template<typename Matrix, typename Vector>
Vector&
SimpleMatrixProvider<Matrix, Vector>::
getVector(MatrixSpecificationsProvider const& msp, std::size_t& id)
{
    auto mat_spec = msp.getMatrixSpecifications();
    return *getVector_<true>(id, mat_spec).first;
    // TODO assert that the returned object always is of the right size
}

template<typename Matrix, typename Vector>
Vector&
SimpleMatrixProvider<Matrix, Vector>::
getVector(Vector const& x)
{
    std::size_t id;
    auto const& res = getVector_<false>(id, x);
    if (!res.second) // no new object has been created
        BLAS::copy(x, *res.first);
    return *res.first;
}

template<typename Matrix, typename Vector>
Vector&
SimpleMatrixProvider<Matrix, Vector>::
getVector(Vector const& x, std::size_t& id)
{
    auto const& res = getVector_<false>(id, x);
    if (!res.second) // no new object has been created
        BLAS::copy(x, *res.first);
    return *res.first;
}

template<typename Matrix, typename Vector>
void
SimpleMatrixProvider<Matrix, Vector>::
releaseVector(Vector const& x)
{
    auto it = _used_vectors.find(const_cast<Vector*>(&x));
    if (it == _used_vectors.end()) {
        ERR("A vector with the id XXX has not been found. Cannot release it. Aborting.");
        // TODO message
        std::abort();
    } else {
        ::detail::transfer(_used_vectors, _unused_vectors, it);
    }
}

template<typename Matrix, typename Vector>
SimpleMatrixProvider<Matrix, Vector>::
~SimpleMatrixProvider()
{
    if ((!_used_matrices.empty()) || (!_used_vectors.empty())) {
        WARN("There are still some matrices and vectors in use."
             " This might be an indicator of a possible waste of memory.");
    }

    for (auto& id_ptr : _unused_matrices)
        delete id_ptr.second;

    for (auto& ptr_id : _used_matrices)
        delete ptr_id.first;

    for (auto& id_ptr : _unused_vectors)
        delete id_ptr.second;

    for (auto& ptr_id : _used_vectors)
        delete ptr_id.first;
}

} // MathLib
