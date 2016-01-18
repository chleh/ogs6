#pragma once

#include <cassert>

namespace BaseLib
{

enum class StorageOrder { ColumnMajor, RowMajor };


template <typename T, long Rows, long Cols>
class MatrixRef
{
public:
    explicit MatrixRef(T* data, StorageOrder order)
        : _data{data}, _order{order} {}

    MatrixRef(MatrixRef<T, Rows, Cols> const& other) = default;
    MatrixRef(MatrixRef<T, Rows, Cols> && other) = default;

    /* operator= might be confusing if used in actual code */
    MatrixRef<T, Rows, Cols> operator=(MatrixRef<T, Rows, Cols> const& other) = delete;
    MatrixRef<T, Rows, Cols> operator=(MatrixRef<T, Rows, Cols> && other) = delete;

    T& operator()(long row, long col)
    {
        assert(row >= 0 && row < Rows && col >= 0 && col < Cols);

        switch (_order)
        {
        case StorageOrder::ColumnMajor:
            return _data[col*Rows + row];
        case StorageOrder::RowMajor:
            return _data[row*Cols + col];
        }
        return *_data; // cannot happen
    }
    T const& operator()(long row, long col) const
    {
        assert(row >= 0 && row < Rows && col >= 0 && col < Cols);

        switch (_order)
        {
        case StorageOrder::ColumnMajor:
            return _data[col*Rows + row];
        case StorageOrder::RowMajor:
            return _data[row*Cols + col];
        }
        return *_data; // cannot happen
    }

    StorageOrder getStorageOrder() const { return _order; }

    long rows() const { return Rows; }
    long cols() const { return Cols; }

private:
    T* _data;
    StorageOrder _order;
};


} // namespace BaseLib
