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
    : _A(A), _rhs(rhs), _data_pos(data_pos) {}

    ~VectorMatrixAssembler() {}

    void setX(GLOBAL_VECTOR_ const * x, GLOBAL_VECTOR_ const * x_prev_ts)
    {
        assert((!x == !x_prev_ts) && "either no or both inputs have to be set");
        assert((!x) || x->size() == x_prev_ts->size());
        _x = x;
        _x_prev_ts = x_prev_ts;
    }


    void getLocalNodalValues(std::size_t const id,
                             std::vector<double>& localX,
                             std::vector<double>& localX_prev_ts) const
    {
        std::vector<GlobalIndexType> indices;
        getLocalNodalValuesIndices(id, localX, localX_prev_ts, indices);
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
        std::vector<double> localX;
        std::vector<double> localX_pts;
        std::vector<GlobalIndexType> indices;

        if (_x != nullptr)
        {
            getLocalNodalValuesIndices(id, localX, localX_pts, indices);
        }

        LocalToGlobalIndexMap::RowColumnIndices const r_c_indices(
                    indices, indices);

        local_assembler->assemble(localX, localX_pts);
        local_assembler->addToGlobal(_A, _rhs, r_c_indices);
    }

private:
    void getLocalNodalValuesIndices(
            std::size_t const id,
            std::vector<double>& localX,
            std::vector<double>& localX_prev_ts,
            std::vector<GlobalIndexType>& indices
            ) const
    {
        assert(_x != nullptr && _x_prev_ts != nullptr
               && _x->size() == _x_prev_ts->size());
        assert(_data_pos.size() > id);

        indices.clear();

        // Local matrices and vectors will always be ordered by component,
        // no matter what the order of the global matrix is.
        for (unsigned c=0; c<_data_pos.getNumComponents(); ++c)
        {
            auto const& idcs = _data_pos(id, c).rows;
            indices.reserve(indices.size() + idcs.size());
            indices.insert(indices.end(), idcs.begin(), idcs.end());
        }

        localX.reserve(indices.size());
        localX_prev_ts.reserve(indices.size());

        for (auto i : indices)
        {
            localX.emplace_back(_x->get(i));
            localX_pts.emplace_back(_x_prev_ts->get(i));
        }
    }

protected:
    GLOBAL_MATRIX_ &_A;
    GLOBAL_VECTOR_ &_rhs;
    GLOBAL_VECTOR_ const *_x = nullptr;
    GLOBAL_VECTOR_ const *_x_prev_ts = nullptr;
    LocalToGlobalIndexMap const& _data_pos;
};

}   // namespace AssemblerLib

#endif  // ASSEMBLERLIB_VECTORMATRIXASSEMBLER_H_
