#pragma once

#include <memory>

#include "boost/property_tree/ptree.hpp"

#include "OdeSolver.h"

#include "CVodeSolver.h"

namespace MathLib
{


template <unsigned NumEquations>
std::unique_ptr<OdeSolver<NumEquations> >
createOdeSolver(const boost::property_tree::ptree& config);


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

    void setFunction(Function f, JacobianFunction df) override {
        Implementation::setFunction(f, df);
    }

    void setIC(const double t0, const Arr& y0) override {
        Implementation::setIC(t0, y0.data());
    }

    void preSolve() {
        Implementation::preSolve();
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

private:
    /// instances of this class shall only be constructed by
    /// the friend function listed below
    ConcreteOdeSolver(typename Implementation::ConfigTree const& config)
        : Implementation{config}
    {}

    friend std::unique_ptr<OdeSolver<NumEquations> >
    createOdeSolver<NumEquations>(const boost::property_tree::ptree& config);
};


template <unsigned NumEquations>
std::unique_ptr<OdeSolver<NumEquations> >
createOdeSolver(const boost::property_tree::ptree& config)
{
    auto up = std::unique_ptr<OdeSolver<NumEquations> >();
    auto p  = new ConcreteOdeSolver<NumEquations, CVodeSolverInternal>(config);
    up.reset(p);
    return up;
}

} // namespace MathLib
