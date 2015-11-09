/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>
#include <cassert>


namespace MathLib
{

/**
 * Global matrix based on Eigen sparse matrix
 *
 * The matrix will be dynamically allocated during construction.
 */
class CSRMatrix final
{
public:
    using IndexType = int;

    explicit CSRMatrix(IndexType num_cols,
                       std::vector<double>&& values,
                       std::vector<IndexType>&& row_idcs,
                       std::vector<IndexType>&& col_idcs)
        : _num_cols{num_cols}
        , _values{std::move(values)}
        , _row_idcs{std::move(row_idcs)}
        , _col_idcs{std::move(col_idcs)}
    {
        assert(_values.size() == _col_idcs.size());
        assert(_row_idcs.size() >= 1);
        assert(_row_idcs.back() == _values.size());
    }

    /// return the number of rows
    std::size_t getNRows() const { return _row_idcs.size() - 1; }

    /// return the number of columns
    std::size_t getNCols() const { return _num_cols; }

    std::vector<double> const& getValues() const { return _values; }
    std::vector<IndexType> const& getRowIndices() const { return _row_idcs; }
    std::vector<IndexType> const& getColumnIndices() const { return _col_idcs; }


private:
    IndexType _num_cols;
    std::vector<double> _values;
    std::vector<IndexType> _row_idcs;
    std::vector<IndexType> _col_idcs;

};




/*
template <class T_DENSE_MATRIX>
void COOMatrix::add(std::vector<IndexType> const& row_pos,
                      std::vector<IndexType> const& col_pos,
                      const T_DENSE_MATRIX& sub_matrix, double fkt)
{
    auto const n_rows = row_pos.size();
    auto const n_cols = col_pos.size();
    for (auto i = decltype(n_rows){0}; i < n_rows; i++) {
        auto const row = row_pos[i];
        for (auto j = decltype(n_cols){0}; j < n_cols; j++) {
            auto const col = col_pos[j];
            add(row, col, fkt * sub_matrix(i, j));
        }
    }
}
*/

} // end namespace MathLib
