#pragma once

namespace ProcessLib
{

namespace Ode
{

typedef void (* const Function)(const double t, double const*const y, double *const ydot);

class OdeSolver
{
public:
    virtual void init(const unsigned num_equations) = 0;
    virtual void solve(Function f, const double t ) = 0;

    virtual ~OdeSolver() = default;
};

class CVodeSolverImpl;

class CVodeSolver : public OdeSolver
{
public:
    CVodeSolver();
    void init(const unsigned num_equations) override;

    void setTolerance(const double* abstol, const double reltol);
    void setTolerance(const double abstol, const double reltol);

    void setIC(const double t0, double const*const y0);

    void solve(Function f, const double t) override;

    ~CVodeSolver();

private:
    CVodeSolverImpl* _impl;
};


}

}
