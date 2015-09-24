#ifndef MATHLIB_ODESOLVER_H
#define MATHLIB_ODESOLVER_H

#include <array>

#include "declarations.h"

namespace MathLib
{

/**
 * ODE solver Interface
 */
template<unsigned NumEquations, typename... FunctionArguments>
class OdeSolver
{
public:
    using Arr = std::array<double, NumEquations>;
    using Function = MathLib::Function<FunctionArguments...>;
    using JacobianFunction = MathLib::JacobianFunction<FunctionArguments...>;

    virtual void init() = 0;

    virtual void setTolerance(const Arr& abstol, const double reltol) = 0;
    virtual void setTolerance(const double abstol, const double reltol) = 0;

    virtual void setFunction(Function f, JacobianFunction df,
                             FunctionArguments*... args) = 0;

    virtual void setIC(const double t0, const Arr& y0) = 0;

    virtual void preSolve() = 0;
    virtual void solve(const double t) = 0;

    virtual unsigned getNumEquations() const { return NumEquations; }

    virtual double const* getSolution() const = 0;
    virtual double getTime() const = 0;

    virtual ~OdeSolver() = default;
};

} // namespace MathLib

#endif // MATHLIB_ODESOLVER_H
