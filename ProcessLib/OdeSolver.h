#pragma once

#include <array>

namespace ProcessLib
{

namespace Ode
{

// maybe use Eigen::Map here
typedef void (* const Function)(const double t, double const*const y, double *const ydot);

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

    virtual ~OdeSolver() = default;
};


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

    ~CVodeSolverInternal();
private:
    CVodeSolverImpl* _impl; ///< pimpl idiom hides implementation
};


/**
 * ODE solver with a bounds-safe interface
 */
template<unsigned NumEquations>
class CVodeSolver
        : public OdeSolver<NumEquations>,
        private CVodeSolverInternal
{
public:
    using Arr = typename OdeSolver<NumEquations>::Arr;

    void init() override {
        CVodeSolverInternal::init(NumEquations);
    }

    void setTolerance(const Arr& abstol, const double reltol) override {
        CVodeSolverInternal::setTolerance(abstol.data(), reltol);
    }

    void setTolerance(const double abstol, const double reltol) override {
        CVodeSolverInternal::setTolerance(abstol, reltol);
    }

    void setIC(const double t0, const Arr& y0) override {
        CVodeSolverInternal::setIC(t0, y0.data());
    }

    void solve(Function f, const double t) override {
        CVodeSolverInternal::solve(f, t);
    }
};


}

}

#include "OdeSolver-impl.h"
