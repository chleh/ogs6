extern "C"
{
#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <cvode/cvode_dense.h>       /* prototype for CVDense */
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h> /* definition of type realtype */
}

#include <cassert>

#include "logog/include/logog.hpp"

#include "CVodeSolver.h"

namespace
{

struct UserData
{
    MathLib::Function f = nullptr;
    MathLib::JacobianFunction df = nullptr;
};

}

namespace MathLib
{

class CVodeSolverImpl
{
    static_assert(std::is_same<realtype, double>::value, "cvode's realtype is not the same as double");
public:
    CVodeSolverImpl() = default;
    void init(const unsigned num_equations);

    void setTolerance(const double* abstol, const double reltol);
    void setTolerance(const double abstol, const double reltol);

    void setIC(const double t0, double const*const y0);

    void setFunction(Function f, JacobianFunction df);

    void solve(const double t_end);

    double const* getSolution() const { return NV_DATA_S(_y); }
    double getTime() const { return _t; }

    ~CVodeSolverImpl();

private:
    N_Vector _y = nullptr;
    N_Vector _ydot = nullptr;

    realtype _t;

    N_Vector _abstol = nullptr;
    realtype _reltol;

    unsigned _num_equations;
    void*    _cvode_mem;

    UserData _data;
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

void CVodeSolverImpl::setFunction(Function f, JacobianFunction df)
{
    _data.f = f;
    _data.df = df;
}

void CVodeSolverImpl::setIC(const double t0, double const*const y0)
{
    for (unsigned i=0; i<_num_equations; ++i)
    {
        NV_Ith_S(_y, i) = y0[i];
    }

    _t = t0;
}

void CVodeSolverImpl::solve(const double t_end)
{
	assert(_data.f != nullptr && "ode function y'=f(t,y) was not provided");

	auto f_wrapped
			= [](const realtype t, const N_Vector y, N_Vector ydot, void* data)
	{
		((UserData*) data)->f(t, NV_DATA_S(y), NV_DATA_S(ydot));
		return 0;
	};


	int flag = CVodeInit(_cvode_mem, f_wrapped, _t, _y);
	// if (check_flag(&flag, "CVodeInit", 1)) return(1);

	flag = CVodeSetUserData(_cvode_mem, (void*) &_data);

	/* Call CVodeSVtolerances to specify the scalar relative tolerance
	 * and vector absolute tolerances */
	flag = CVodeSVtolerances(_cvode_mem, _reltol, _abstol);
	// if (check_flag(&flag, "CVodeSVtolerances", 1)) return(1);

	/* Call CVDense to specify the CVDENSE dense linear solver */
	flag = CVDense(_cvode_mem, _num_equations);
	// if (check_flag(&flag, "CVDense", 1)) return(1);

	if (_data.df)
	{
		auto df_wrapped
				= [](const long /*N*/, const realtype t,
				  const N_Vector y, const N_Vector ydot,
				  const DlsMat jac, void* data,
				  N_Vector /*tmp1*/, N_Vector /*tmp2*/, N_Vector /*tmp3*/
				  ) -> int
		{
			// Caution: by calling the DENSE_COL() macro we assume that matrices
			//          are stored contiguous in memory!
			((UserData*) data)->df(t, NV_DATA_S(y), NV_DATA_S(ydot),
								   DENSE_COL(jac, 0), StorageOrder::ColumnMajor);
			return 1;
		};

		flag = CVDlsSetDenseJacFn(_cvode_mem, df_wrapped);
		// if (check_flag(&flag, "CVDlsSetDenseJacFn", 1)) return(1);
	}

	realtype t_reached;
	flag = CVode(_cvode_mem, t_end, _y, &t_reached, CV_NORMAL);
	_t = t_reached;
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

void CVodeSolverInternal::setFunction(Function f, JacobianFunction df)
{
    _impl->setFunction(f, df);
}

void CVodeSolverInternal::setIC(const double t0, double const*const y0)
{
    _impl->setIC(t0, y0);
}

void CVodeSolverInternal::solve(const double t_end)
{
	_impl->solve(t_end);
}

double const* CVodeSolverInternal::getSolution() const
{
	return _impl->getSolution();
}

double CVodeSolverInternal::getTime() const
{
	return _impl->getTime();
}

CVodeSolverInternal::~CVodeSolverInternal()
{
    delete _impl;
}

} // namespace MathLib

