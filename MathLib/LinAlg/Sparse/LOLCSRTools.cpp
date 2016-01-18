#include "LOLCSRTools.h"



namespace MathLib
{

CSRMatrix
toCSR(LOLMatrix const& mat)
{
    using IT = LOLMatrix::IndexType;

    auto const& data = mat.getData();

    IT nnz = 0;
    for (auto const& row : data) nnz += row.size();

    std::vector<double> values;
    std::vector<IT> row_idcs;
    std::vector<IT> col_idcs;

    values.reserve(nnz);
    col_idcs.reserve(nnz);
    row_idcs.reserve(data.size()+1);

    IT row_idx = 0;
    for (auto const& row : data) {
        row_idcs.push_back(row_idx);
        row_idx += row.size();

        for (auto const& e : row) {
            values.push_back(e.value);
            col_idcs.push_back(e.col_idx);
        }
    }
    row_idcs.push_back(row_idx);

    return CSRMatrix(mat.getNCols(), std::move(values),
                     std::move(row_idcs), std::move(col_idcs));
}


LOLMatrix
toLOL(CSRMatrix const& mat)
{
    using Row = LOLMatrix::Row;
    using IT = LOLMatrix::IndexType;

    const IT nrows = mat.getNRows();
    std::vector<Row> data(nrows);

    auto const& row_idcs = mat.getRowIndices();
    auto const& col_idcs = mat.getColumnIndices();
    auto const& values = mat.getValues();

    for (IT r=0; r<nrows; ++r) {
        const IT lb = row_idcs[r];
        const IT ub = row_idcs[r+1];
        if (lb == ub) continue;
        assert(lb < ub);

        auto& row = data[r];
        row.reserve(ub-lb);

        for (IT i=lb; i<ub; ++i) {
            row.emplace_back(col_idcs[i], values[i]);
        }
    }

    return LOLMatrix(mat.getNCols(), std::move(data));
}

}

