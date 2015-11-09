/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <cassert>
#include <vector>
#include <algorithm>

#include "MathLib/LinAlg/RowColumnIndices.h"

namespace MathLib
{

/**
 * Global matrix based on Eigen sparse matrix
 *
 * The matrix will be dynamically allocated during construction.
 */
class LOLMatrix final
{
public:
    using IndexType = int; // TODO: use mkl int or so

    struct Entry {
        Entry(IndexType i, double v) : col_idx(i), value(v) {}

        IndexType col_idx;
        double value;

        friend bool operator< (Entry const& self, Entry const& other) {
            return self.col_idx < other.col_idx;
        }
        friend bool operator< (Entry const& e, IndexType idx) {
            return e.col_idx < idx;
        }
        friend bool operator< (IndexType idx, Entry const& e) {
            return idx < e.col_idx;
        }
        friend bool operator== (Entry const& a, Entry const& b) {
            return a.col_idx == b.col_idx
                    && a.value == b.value;
        }
    };

    using Row = std::vector<Entry>;

    /**
     * constructor
     * @param n the number of rows (that is equal to the number of columns)
     * @param n_nonzero_columns the number of non-zero columns used for preallocation
     */
    explicit LOLMatrix(IndexType num_rows, IndexType num_cols)
        : _mat{num_rows}, _num_cols{num_cols}
    {}
    explicit LOLMatrix(IndexType num_rows)
        : _mat{num_rows}, _num_cols{num_rows}
    {}


    explicit LOLMatrix(IndexType num_cols, std::vector<Row>&& data)
        : _mat{std::move(data)}, _num_cols{num_cols}
    {
#ifndef NDEBUG
        for (auto const& row : _mat) {
            for (auto const& e : row) {
                assert(e.col_idx < _num_cols);
            }
        }
#endif
    }

    /// return the number of rows
    IndexType getNRows() const { return _mat.size(); }

    /// return the number of columns
    IndexType getNCols() const { return _num_cols; }

    /// reset data entries to zero.
    void setZero()
    {
        for (auto& row : _mat) {
            for (auto& e : row) {
                e.value = 0.0;
            }
        }
    }

    /// set a value to the given entry. If the entry doesn't exist, this class
    /// dynamically allocates it.
    void setValue(IndexType row_idx, IndexType col_idx, double value)
    {
        assert(row_idx >= 0 && row_idx < _mat.size());

        auto& row = _mat[row_idx];

        auto end = row.end();
        auto pos = std::lower_bound(row.begin(), end, col_idx);

        if (pos == end) {
            row.emplace_back(col_idx, value);
        } else if (pos->col_idx == col_idx) {
            pos->value = value;
        } else {
            row.emplace(pos, col_idx, value);
        }

        assert_sorted(row);
    }

    /// add a value to the given entry. If the entry doesn't exist, the value is
    /// inserted.
    int add(IndexType row_idx, IndexType col_idx, double value)
    {
        assert(row_idx >= 0 && row_idx < _mat.size());

        auto& row = _mat[row_idx];

        auto end = row.end();
        auto pos = std::lower_bound(row.begin(), end, col_idx);

        if (pos == end) {
            row.emplace_back(col_idx, value);
        } else if (pos->col_idx == col_idx) {
            pos->value += value;
        } else {
            row.emplace(pos, col_idx, value);
        }

        assert_sorted(row);
    }

    /// Add sub-matrix at positions \c row_pos and same column positions as the
    /// given row positions. If the entry doesn't exist, the value is inserted.
    template<class T_DENSE_MATRIX>
    void add(std::vector<IndexType> const& row_pos,
            const T_DENSE_MATRIX &sub_matrix,
            double fkt = 1.0)
    {
        this->add(row_pos, row_pos, sub_matrix, fkt);
    }

    /// Add sub-matrix at positions given by \c indices. If the entry doesn't exist,
    /// this class inserts the value.
    template<class T_DENSE_MATRIX>
    void add(RowColumnIndices<IndexType> const& indices,
            const T_DENSE_MATRIX &sub_matrix,
            double fkt = 1.0)
    {
        this->add(indices.rows, indices.columns, sub_matrix, fkt);
    }

    /// Add sub-matrix at positions \c row_pos and \c col_pos. If the entries doesn't
    /// exist in the matrix, the values are inserted.
    /// @param row_pos     a vector of row position indices. The vector size should
    ///                    equal to the number of rows in the given sub-matrix.
    /// @param col_pos     a vector of column position indices. The vector size should
    ///                    equal to the number of columns in the given sub-matrix.
    /// @param sub_matrix  a sub-matrix to be added
    /// @param fkt         a scaling factor applied to all entries in the sub-matrix
    template <class T_DENSE_MATRIX>
    void add(std::vector<IndexType> const& row_pos,
            std::vector<IndexType> const& col_pos, const T_DENSE_MATRIX &sub_matrix,
            double fkt = 1.0);

    /// get value. This function returns zero if the element doesn't exist.
    double get(IndexType row_idx, IndexType col_idx) const
    {
        assert(row_idx >= 0 && row_idx < _mat.size());

        auto& row = _mat[row_idx];

        auto end = row.end();
        auto pos = std::lower_bound(row.begin(), end, col_idx);

        if (pos == end || pos->col_idx != col_idx) {
            return 0.0;
        } else if (pos->col_idx == col_idx) {
            return pos->value;
        }
    }

    /// get value. This function returns zero if the element doesn't exist.
    double operator() (IndexType row, IndexType col) const
    {
        return get(row, col);
    }

#if 0
    /// get a maximum value in diagonal entries
    double getMaxDiagCoeff() const
    {
        return _mat.diagonal().maxCoeff();
    }

    /// y = mat * x
    void multiply(const EigenVector &x, EigenVector &y) const
    {
        y.getRawVector() = _mat * x.getRawVector();
    }
#endif

    /// return always true, i.e. the matrix is always ready for use
    bool isAssembled() const { return true; }

    std::vector<Row> const& getData() const { return _mat; }

    friend bool operator==(LOLMatrix const& a, LOLMatrix const& b)
    {
        return a._num_cols == b._num_cols
                && a._mat == b._mat;
    }

    void write(std::string const& /*path*/)
    {
        // TODO implement
    }

protected:
    std::vector<Row> _mat;
    const IndexType _num_cols;


    void assert_sorted(Row const& row)
    {
#ifndef NDEBUG
        for (IndexType i=0; i<row.size()-1; ++i)
            assert(row[i] < row[i+1]);
#endif
    }
};



template <class T_DENSE_MATRIX>
void LOLMatrix::add(std::vector<IndexType> const& row_pos,
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

} // end namespace MathLib



