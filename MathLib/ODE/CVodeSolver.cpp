extern "C"
{
#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <cvode/cvode_dense.h>       /* prototype for CVDense */
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h> /* definition of type realtype */
}

#include "logog/include/logog.hpp"

#include "CVodeSolver.h"

namespace MathLib
{

class CVodeSolverImpl
{
public:
    CVodeSolverImpl() = default;
    void init(const unsigned num_equations);

    void setTolerance(const double* abstol, const double reltol);
    void setTolerance(const double abstol, const double reltol);

    void setIC(const double t0, double const*const y0);

    void solve(Function f, const double t);

    double const* getSolution() const { return NV_DATA_S(_y); }

    ~CVodeSolverImpl();

private:
    N_Vector _y = nullptr;
    N_Vector _ydot = nullptr;

    realtype _t;

    N_Vector _abstol = nullptr;
    realtype _reltol;

    unsigned _num_equations;
    void*    _cvode_mem;

};


void CVodeSolverImpl::init(const unsigned num_equations)
{
    _y      = N_VNew_Serial(num_equations);
    _ydot   = N_VNew_Serial(num_equations);
    _abstol = N_VNew_Serial(num_equations);
    _num_equations = num_equations;

    _cvode_mem = CVodeCreate(CV_ADAMS, CV_FUNCTIONAL);
}

void CVodeSolverImpl::setTolerance(const double *abstol, const double reltol)
{
    for (unsigned i=0; i<_num_equations; ++i)
    {
        NV_Ith_S(_abstol, i) = abstol[i];
    }

    _reltol = reltol;
}

void CVodeSolverImpl::setTolerance(const double abstol, const double reltol)
{
    for (unsigned i=0; i<_num_equations; ++i)
    {
        NV_Ith_S(_abstol, i) = abstol;
    }

    _reltol = reltol;
}

void CVodeSolverImpl::setIC(const double t0, double const*const y0)
{
    for (unsigned i=0; i<_num_equations; ++i)
    {
        NV_Ith_S(_y, i) = y0[i];
    }

    _t = t0;
}

void CVodeSolverImpl::solve(Function f, const double t)
{
	auto f_wrapped
			= [](const realtype t, const N_Vector y, N_Vector ydot, void* f)
	{
		((Function) f) (t, NV_DATA_S(y), NV_DATA_S(ydot));
		return 0;
	};


	int flag = CVodeInit(_cvode_mem, f_wrapped, _t, _y);
	// if (check_flag(&flag, "CVodeInit", 1)) return(1);

	flag = CVodeSetUserData(_cvode_mem, (void*) f);

	/* Call CVodeSVtolerances to specify the scalar relative tolerance
	 * and vector absolute tolerances */
	flag = CVodeSVtolerances(_cvode_mem, _reltol, _abstol);
	// if (check_flag(&flag, "CVodeSVtolerances", 1)) return(1);

	/* Call CVDense to specify the CVDENSE dense linear solver */
	flag = CVDense(_cvode_mem, _num_equations);
	// if (check_flag(&flag, "CVDense", 1)) return(1);

	realtype t_end;
	flag = CVode(_cvode_mem, t, _y, &t_end, CV_NORMAL);
	// std::cout << "result at time " << t << " is " << NV_Ith_S(y,0) << std::endl;
	if (flag != CV_SUCCESS) {
		// std::cerr << "ERROR at " << __FUNCTION__ << ":" << __LINE__ << std::endl;
	}
}

CVodeSolverImpl::~CVodeSolverImpl()
{
    if (_y) {
        N_VDestroy_Serial(_y);
        N_VDestroy_Serial(_ydot);
        N_VDestroy_Serial(_abstol);
    }

    if (_cvode_mem)
    {
        CVodeFree(&_cvode_mem);
    }
}




CVodeSolverInternal::CVodeSolverInternal()
    : _impl(new CVodeSolverImpl)
{}

void CVodeSolverInternal::init(const unsigned num_equations)
{
    _impl->init(num_equations);
}

void CVodeSolverInternal::setTolerance(const double *abstol, const double reltol)
{
    _impl->setTolerance(abstol, reltol);
}

void CVodeSolverInternal::setTolerance(const double abstol, const double reltol)
{
    _impl->setTolerance(abstol, reltol);
}

void CVodeSolverInternal::setIC(const double t0, double const*const y0)
{
    _impl->setIC(t0, y0);
}

void CVodeSolverInternal::solve(Function f, const double t)
{
	_impl->solve(f, t);
}

double const* CVodeSolverInternal::getSolution() const
{
	return _impl->getSolution();
}

CVodeSolverInternal::~CVodeSolverInternal()
{
    delete _impl;
}

} // namespace MathLib

