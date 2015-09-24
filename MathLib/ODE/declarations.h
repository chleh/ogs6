#pragma once

#include<cassert>

#include "BaseLib/ArrayRef.h"

namespace MathLib
{

enum class StorageOrder { ColumnMajor, RowMajor };


template<unsigned N, typename... FunctionArguments>
using Function = bool (*)(const double t,
                          BaseLib::ArrayRef<const double, N> y,
                          BaseLib::ArrayRef<double, N> ydot,
                          FunctionArguments&... arg);

template<unsigned N, typename... FunctionArguments>
using JacobianFunction = bool (*)(const double t,
                                  BaseLib::ArrayRef<const double, N> y,
                                  BaseLib::ArrayRef<const double, N> ydot,
                                  double *const jac, // TODO: write matrix wrapper class
                                  StorageOrder order,
                                  FunctionArguments&... arg);


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



// This is an internal detail
class FunctionHandles
{
public:
    virtual bool call(const double t, double const*const y, double *const ydot
                      ) = 0;
    virtual bool callJacobian(
            const double t,
            double const*const y,
            double const*const ydot,
            double *const jac,
            StorageOrder order
            ) = 0;

    virtual bool hasJacobian() const = 0;

    virtual unsigned getNumEquations() const = 0;

    virtual ~FunctionHandles() = default;
};



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
