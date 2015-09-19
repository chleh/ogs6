#pragma once

#include <memory>

#include "OdeSolver.h"

#include "CVodeSolver.h"

namespace MathLib
{

/**
 * ODE solver with a bounds-safe interface
 */
template<unsigned NumEquations, typename Implementation>
class ConcreteOdeSolver
        : public OdeSolver<NumEquations>,
        private Implementation
{
public:
    using Arr = typename OdeSolver<NumEquations>::Arr;

    void init() override {
        Implementation::init(NumEquations);
    }

    void setTolerance(const Arr& abstol, const double reltol) override {
        Implementation::setTolerance(abstol.data(), reltol);
    }

    void setTolerance(const double abstol, const double reltol) override {
        Implementation::setTolerance(abstol, reltol);
    }

    void setFunction(Function f) override {
        Implementation::setFunction(f);
    }

    void setIC(const double t0, const Arr& y0) override {
        Implementation::setIC(t0, y0.data());
    }

    void solve(const double t) override {
        Implementation::solve(t);
    }

    double const* getSolution() const override {
        return Implementation::getSolution();
    }

    double getTime() const override {
        return Implementation::getTime();
    }
};


template <unsigned NumEquations>
std::unique_ptr<OdeSolver<NumEquations> > createOdeSolver()
{
    auto up = std::unique_ptr<OdeSolver<NumEquations> >();
    auto p  = new ConcreteOdeSolver<NumEquations, CVodeSolverInternal>;
    up.reset(p);
    return up;
}

} // namespace MathLib
