// TODO doc

#include <cassert>
#include <logog/include/logog.hpp>

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
    // return res.first->second;
}

} // detail


namespace MathLib
{

template<typename Matrix, typename Vector>
Matrix&
SimpleMatrixProvider<Matrix, Vector>::
getMatrix()
{
    // not found, so create a new one
    auto const id = _next_id++;
    auto res = _used_matrices.emplace(new Matrix(0, 0), id);
    assert(res.second && "Emplacement failed.");
    return *res.first->first;
}

template<typename Matrix, typename Vector>
Matrix&
SimpleMatrixProvider<Matrix, Vector>::
getMatrix(MatrixSpecificationsProvider const& msp, std::size_t& id)
{
    auto it = _unused_matrices.find(id);
    if (it == _unused_matrices.end())
    {
        auto mat_spec = msp.getMatrixSpecifications();
        // not found, so create a new one
        id = _next_id++;
        auto res = _used_matrices.emplace(
            new Matrix{mat_spec.nrows, mat_spec.ncols}, id);
        assert(res.second && "Emplacement failed.");
        return *res.first->first;
    }
    else { // unused matrix found
        return *::detail::transfer(_unused_matrices, _used_matrices, it);
    }
}

template<typename Matrix, typename Vector>
Matrix&
SimpleMatrixProvider<Matrix, Vector>::
getMatrix(MatrixSpecificationsProvider const& msp, std::size_t& id, Matrix const& A)
{
    (void) msp; // TODO
    auto it = _unused_matrices.find(id);
    if (it == _unused_matrices.end())
    {
        // not found, so create a new one
        id = _next_id++;
        auto res = _used_matrices.emplace(new Matrix{A}, id);
        assert(res.second && "Emplacement failed.");
        return *res.first->first;
    }
    else { // unused matrix found
        return *::detail::transfer(_unused_matrices, _used_matrices, it);
    }
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
Vector&
SimpleMatrixProvider<Matrix, Vector>::
getVector()
{
    // TODO create an empty zero-size vector
    std::size_t const id = _next_id++;
    auto res = _used_vectors.emplace(new Vector, id);
    assert(res.second && "Emplacement failed.");
    return *res.first->first;
}

template<typename Matrix, typename Vector>
Vector&
SimpleMatrixProvider<Matrix, Vector>::
getVector(Vector const& x)
{
    // TODO create an empty zero-size vector
    std::size_t const id = _next_id++;
    auto res = _used_vectors.emplace(new Vector{x}, id);
    assert(res.second && "Emplacement failed.");
    return *res.first->first;
}

template<typename Matrix, typename Vector>
Vector&
SimpleMatrixProvider<Matrix, Vector>::
getVector(std::size_t& id)
{
    auto it = _unused_vectors.find(id);
    if (it == _unused_vectors.end())
    {
        // not found, so create a new one
        id = _next_id++;
        auto res = _used_vectors.emplace(new Vector, id);
        assert(res.second && "Emplacement failed.");
        return *res.first->first;
    }
    else { // unused vector found
        return *::detail::transfer(_unused_vectors, _used_vectors, it);
    }
}

template<typename Matrix, typename Vector>
Vector&
SimpleMatrixProvider<Matrix, Vector>::
getVector(MatrixSpecificationsProvider const& msp, std::size_t& id)
{
    auto it = _unused_vectors.find(id);
    if (it == _unused_vectors.end())
    {
        // not found, so create a new one
        id = _next_id++;
        auto res = _used_vectors.emplace(
            new Vector{msp.getMatrixSpecifications().nrows}, id);
        assert(res.second && "Emplacement failed.");
        return *res.first->first;
    }
    else { // unused vector found
        return *::detail::transfer(_unused_vectors, _used_vectors, it);
    }
}

template<typename Matrix, typename Vector>
Vector&
SimpleMatrixProvider<Matrix, Vector>::
getVector(MatrixSpecificationsProvider const& msp, std::size_t& id, Vector const& x)
{
    (void) msp; // TODO
    auto it = _unused_vectors.find(id);
    if (it == _unused_vectors.end())
    {
        // not found, so create a new one
        id = _next_id++;
        auto res = _used_vectors.emplace(new Vector{x}, id);
        assert(res.second && "Emplacement failed.");
        return *res.first->first;
    }
    else { // unused vector found
        return *::detail::transfer(_unused_vectors, _used_vectors, it);
    }
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
