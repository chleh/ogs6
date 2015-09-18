#ifndef MATHLIB_CVODESOLVER_H
#define MATHLIB_CVODESOLVER_H

#include "declarations.h"

namespace MathLib
{

class CVodeSolverImpl;


/**
 * ODE solver, general, pointer based implementation. No implicit bounds checking
 *
 * For internal use only.
 */
class CVodeSolverInternal
{
protected:
    CVodeSolverInternal();
    void init(const unsigned num_equations);

    void setTolerance(double const*const abstol, const double reltol);
    void setTolerance(const double abstol, const double reltol);

    void setIC(const double t0, double const*const y0);

    void solve(Function f, const double t);

    double const* getSolution() const;

    ~CVodeSolverInternal();
private:
    CVodeSolverImpl* _impl; ///< pimpl idiom hides implementation
};

} // namespace MathLib

#endif // MATHLIB_CVODESOLVER_H
