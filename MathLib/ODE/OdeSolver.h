#ifndef MATHLIB_ODESOLVER_H
#define MATHLIB_ODESOLVER_H

#include <array>

#include "declarations.h"

namespace MathLib
{

/**
 * ODE solver Interface
 */
template<unsigned NumEquations>
class OdeSolver
{
public:
    using Arr = std::array<double, NumEquations>;

    virtual void init() = 0;

    virtual void setTolerance(const Arr& abstol, const double reltol) = 0;
    virtual void setTolerance(const double abstol, const double reltol) = 0;

    virtual void setIC(const double t0, const Arr& y0) = 0;

    virtual void solve(Function f, const double t ) = 0;

    virtual unsigned getNumEquations() const { return NumEquations; }

    virtual double const* getSolution() const = 0;

    virtual ~OdeSolver() = default;
};

} // namespace MathLib

#endif // MATHLIB_ODESOLVER_H
