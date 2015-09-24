#pragma once

#include <memory>

#include "boost/property_tree/ptree.hpp"

#include "OdeSolver.h"

#include "CVodeSolver.h"

namespace MathLib
{

namespace detail
{


template<typename... FunctionArguments>
struct Handles;

template<typename FunctionArgument>
struct Handles<FunctionArgument>
        : public MathLib::FunctionHandles
{
    using Function = MathLib::Function<FunctionArgument>;
    using JacobianFunction = MathLib::JacobianFunction<FunctionArgument>;

    bool call(const double t, const double * const y, double * const ydot) override
    {
        if (f) return f(t, y, ydot, *_data);
        return true;
    }

    bool callJacobian(const double t, const double * const y, const double * const ydot,
                      double * const jac, StorageOrder order) override
    {
        if (df) return df(t, y, ydot, jac, order, *_data);
        return true;
    }

    bool hasJacobian() { return df != nullptr; }

    void setArguments(FunctionArgument* arg) { _data = arg; }

    Function f = nullptr;
    JacobianFunction df = nullptr;

private:
    FunctionArgument* _data = nullptr;
};

template<>
struct Handles<> : public MathLib::FunctionHandles
{
    using Function = MathLib::Function<>;
    using JacobianFunction = MathLib::JacobianFunction<>;

    bool call(const double t, const double * const y, double * const ydot) override
    {
        if (f) return f(t, y, ydot);
        return true;
    }

    bool callJacobian(const double t, const double * const y, const double * const ydot,
                      double * const jac, StorageOrder order) override
    {
        if (df) return df(t, y, ydot, jac, order);
        return true;
    }

    bool hasJacobian() { return df != nullptr; }

    void setArguments() const {}

    Function f = nullptr;
    JacobianFunction df = nullptr;
};


} // namespace detail


template <unsigned NumEquations, typename... FunctionArguments>
std::unique_ptr<OdeSolver<NumEquations, FunctionArguments...> >
createOdeSolver(const boost::property_tree::ptree& config);


/**
 * ODE solver with a bounds-safe interface
 */
template<unsigned NumEquations, typename Implementation, typename... FunctionArguments>
class ConcreteOdeSolver
        : public OdeSolver<NumEquations, FunctionArguments...>,
        private Implementation
{
public:
    using Interface = OdeSolver<NumEquations, FunctionArguments...>;
    using Arr = typename Interface::Arr;
    using Function = typename Interface::Function;
    using JacobianFunction = typename Interface::JacobianFunction;

    void init() override {
        Implementation::init(NumEquations);
        Implementation::setFunction(&_handles);
    }

    void setTolerance(const Arr& abstol, const double reltol) override {
        Implementation::setTolerance(abstol.data(), reltol);
    }

    void setTolerance(const double abstol, const double reltol) override {
        Implementation::setTolerance(abstol, reltol);
    }

    void setFunction(Function f, JacobianFunction df, FunctionArguments*... args) override {
        _handles.f = f; _handles.df = df;
        _handles.setArguments(args...);
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

    detail::Handles<FunctionArguments...> _handles;

    friend std::unique_ptr<OdeSolver<NumEquations, FunctionArguments...> >
    createOdeSolver<NumEquations, FunctionArguments...>
    (const boost::property_tree::ptree& config);
};


template <unsigned NumEquations, typename... FunctionArguments>
std::unique_ptr<OdeSolver<NumEquations, FunctionArguments...> >
createOdeSolver(const boost::property_tree::ptree& config)
{
    auto up = std::unique_ptr<OdeSolver<NumEquations, FunctionArguments...> >();
    auto p  = new ConcreteOdeSolver<NumEquations, CVodeSolverInternal, FunctionArguments...>
            (config);
    up.reset(p);
    return up;
}

} // namespace MathLib
