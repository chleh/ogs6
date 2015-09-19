#pragma once

#include<cassert>

namespace MathLib
{

enum class StorageOrder { ColumnMajor, RowMajor };

// maybe use Eigen::Map here
// and use std::function
typedef void (*Function)(const double t, double const*const y, double *const ydot);

typedef void (*JacobianFunction)(const double t,
                                 double const*const y,
                                 double const*const ydot,
                                 double *const jac,
                                 StorageOrder order);

inline void
setMatrixValue(double* matrix,
               const unsigned num_rows, const unsigned num_columns,
               const StorageOrder order,
               const unsigned row, const unsigned column, const double value)
{
    assert(row < num_rows && column < num_columns);

    switch (order)
    {
    case StorageOrder::ColumnMajor:
        matrix[column*num_rows + row] = value;
        break;
    case StorageOrder::RowMajor:
        matrix[row*num_columns + column] = value;
        break;
    }
}

/*
// TODO: values array always row majow?
inline void
setMatrixValues(double* matrix,
                const unsigned num_rows, const unsigned num_columns,
                const StorageOrder order,
                const unsigned row, const unsigned column, const double* values)
{
    switch (order)
    {
    case StorageOrder::ColumnMajor:
        for (unsigned c=0; c<num_columns; ++c) {
            for (unsigned r=0; r<num_rows; ++r) {
                matrix[c*num_rows+r] =
            }
        }
        matrix[column*num_rows + row] = value;
        break;
    case StorageOrder::RowMajor:
        matrix[row*num_columns + column] = value;
        break;
    }
}
*/

}
