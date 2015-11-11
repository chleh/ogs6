/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ASSEMBLERLIB_VECTORMATRIXASSEMBLER_H_
#define ASSEMBLERLIB_VECTORMATRIXASSEMBLER_H_

#include "LocalToGlobalIndexMap.h"

namespace AssemblerLib
{

/// Adds result of local assembler into a global vector and a global matrix.
/// The VectorMatrixAssembler executes the local assembler for a given mesh item
/// and adds the local vector and matrix entries into the global vector and
/// the global matrix. The indices in global objects are provided by
/// the LocalToGlobalIndexMap in the construction.
template<
    typename GLOBAL_MATRIX_,
    typename GLOBAL_VECTOR_>
class VectorMatrixAssembler
{
public:
    typedef GLOBAL_MATRIX_ GLOBAL_MATRIX;
    typedef GLOBAL_VECTOR_ GLOBAL_VECTOR;

public:
    VectorMatrixAssembler(
        GLOBAL_MATRIX_ &A,
        GLOBAL_VECTOR_ &rhs,
        LocalToGlobalIndexMap const& data_pos)
    : _A(A), _rhs(rhs), _data_pos(data_pos),
    _cache{new Cache}
    {}

    ~VectorMatrixAssembler() {}

    void setX(GLOBAL_VECTOR_ const * x, GLOBAL_VECTOR_ const * x_prev_ts)
    {
        assert((!x == !x_prev_ts) && "either no or both inputs have to be set");
        assert((!x) || x->size() == x_prev_ts->size());
        _x = x;
        _x_prev_ts = x_prev_ts;
    }


    void getLocalNodalValues(std::size_t const id,
                             std::vector<double> const*& localX,
                             std::vector<double> const*& localX_prev_ts) const
    {
        getLocalNodalValuesIndices(id);
        localX = &_cache->localX;
        localX_prev_ts = &_cache->localX_prev_ts;
    }


    /// Executes local assembler for the given mesh item and adds the result
    /// into the global matrix and vector.
    /// The positions in the global matrix/vector are taken from
    /// the LocalToGlobalIndexMap provided in the constructor at index \c id.
    /// \attention The index \c id is not necesserily the mesh item's id.
    template <typename LocalAssembler_>
    void operator()(std::size_t const id,
        LocalAssembler_* const local_assembler) const
    {
        getLocalNodalValuesIndices(id);

        LocalToGlobalIndexMap::RowColumnIndices const r_c_indices(
                    _cache->indices, _cache->indices);

        local_assembler->assemble(_cache->localX, _cache->localX_prev_ts);
        local_assembler->addToGlobal(_A, _rhs, r_c_indices);
    }

private:
    void getLocalNodalValuesIndices(std::size_t const id) const
    {

        _cache->indices.clear();

        // Local matrices and vectors will always be ordered by component,
        // no matter what the order of the global matrix is.
        for (unsigned c=0; c<_data_pos.getNumComponents(); ++c)
        {
            auto const& idcs = _data_pos(id, c).rows;
            _cache->indices.reserve(_cache->indices.size() + idcs.size());
            _cache->indices.insert(_cache->indices.end(), idcs.begin(), idcs.end());
        }

        if (_x != nullptr)
        {
            assert(_x != nullptr && _x_prev_ts != nullptr
                   && _x->size() == _x_prev_ts->size());
            assert(_data_pos.size() > id);

            _cache->localX.clear();
            _cache->localX_prev_ts.clear();
            _cache->localX.reserve(_cache->indices.size());
            _cache->localX_prev_ts.reserve(_cache->indices.size());

            for (auto i : _cache->indices)
            {
                _cache->localX.emplace_back(_x->get(i));
                _cache->localX_prev_ts.emplace_back(_x_prev_ts->get(i));
            }
        }
    }

protected:
    GLOBAL_MATRIX_ &_A;
    GLOBAL_VECTOR_ &_rhs;
    GLOBAL_VECTOR_ const *_x = nullptr;
    GLOBAL_VECTOR_ const *_x_prev_ts = nullptr;
    LocalToGlobalIndexMap const& _data_pos;

private:
    struct Cache
    {
        // cache for results
        std::vector<double> localX;
        std::vector<double> localX_prev_ts;
        std::vector<GlobalIndexType> indices;
    };

    std::unique_ptr<Cache> _cache;
};

}   // namespace AssemblerLib

#endif  // ASSEMBLERLIB_VECTORMATRIXASSEMBLER_H_
