// TODO doc

#include <cassert>
#include <logog/include/logog.hpp>

#include "SimpleMatrixProvider.h"

namespace detail
{

template<typename MatVec>
std::unique_ptr<MatVec>&
transfer(std::map<std::size_t, std::unique_ptr<MatVec> >& from,
         std::map<std::size_t, std::unique_ptr<MatVec> >& to,
         typename std::map<std::size_t, std::unique_ptr<MatVec> >::iterator it)
{
    //
    auto res = to.emplace(std::move(*it));
    assert(res.second && "Emplacement failed.");
    from.erase(it);
    //return res.first->second;
    return std::get<1>(*res.first);
}

} // detail


namespace MathLib
{

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
        auto res = _used_matrices.emplace(id,
            std::unique_ptr<Matrix>{new Matrix(mat_spec.nrows, mat_spec.ncols)});
        assert(res.second && "Emplacement failed.");
        return *res.first->second;
    }
    else { // unused matrix found
        return *detail::transfer(_unused_matrices, _used_matrices, it);
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
        auto res = _used_matrices.emplace(id,
            std::unique_ptr<Matrix>{new Matrix(A)});
        assert(res.second && "Emplacement failed.");
        return *res.first->second;
    }
    else { // unused matrix found
        return *detail::transfer(_unused_matrices, _used_matrices, it);
    }
}

template<typename Matrix, typename Vector>
void
SimpleMatrixProvider<Matrix, Vector>::
releaseMatrix(std::size_t const id, Matrix const& /*A*/)
{
    auto it = _used_matrices.find(id);
    if (it == _used_matrices.end()) {
        ERR("A matrix with the id %lu has not been found. Cannot release it. Aborting.");
        std::abort();
    } else {
        detail::transfer(_used_matrices, _unused_matrices, it);
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
        auto res = _used_vectors.emplace(id,
            std::unique_ptr<Vector>{new Vector(msp.getMatrixSpecifications().nrows)});
        assert(res.second && "Emplacement failed.");
        return *res.first->second;
    }
    else { // unused vector found
        return *detail::transfer(_unused_vectors, _used_vectors, it);
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
        auto res = _used_vectors.emplace(id,
            std::unique_ptr<Vector>{new Vector(x)});
        assert(res.second && "Emplacement failed.");
        return *res.first->second;
    }
    else { // unused vector found
        return *detail::transfer(_unused_vectors, _used_vectors, it);
    }
}

template<typename Matrix, typename Vector>
void
SimpleMatrixProvider<Matrix, Vector>::
releaseVector(std::size_t const id, Vector const& x)
{
    auto it = _used_vectors.find(id);
    if (it == _used_vectors.end()) {
        ERR("A vector with the id %lu has not been found. Cannot release it. Aborting.", id);
        std::abort();
    } else {
        detail::transfer(_used_vectors, _unused_vectors, it);
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
}

} // MathLib
