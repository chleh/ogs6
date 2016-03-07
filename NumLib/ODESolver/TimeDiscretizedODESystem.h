/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef NUMLIB_TIMEDISCRETIZEDODESYSTEM_H
#define NUMLIB_TIMEDISCRETIZEDODESYSTEM_H

#include <memory>

#include "MathLib/LinAlg/ApplyKnownSolution.h"
#include "ProcessLib/DirichletBc.h"

#include "ODESystem.h"
#include "NonlinearSystem.h"
#include "TimeDiscretization.h"
#include "MatrixTranslator.h"


namespace NumLib
{

//! \addtogroup ODESolver
//! @{

/*! A NonlinearSystem together with some TimeDiscretization scheme.
 *
 * This is the interface of an ODE towards the TimeLoop.
 * This interface is abstract, it represents any type of first-order ODE.
 *
 * \tparam Matrix the type of matrices occuring in the linearization of the ODE.
 * \tparam Vector the type of the solution vector of the ODE.
 * \tparam NLTag  a tag indicating the method used for resolving nonlinearities.
 */
template<typename Matrix, typename Vector, NonlinearSolverTag NLTag>
class TimeDiscretizedODESystemBase
        : public NonlinearSystem<Matrix, Vector, NLTag>
        , public InternalMatrixStorage
{
public:
    //! Exposes the used time discretization scheme.
    virtual TimeDiscretization<Vector>& getTimeDiscretization() = 0;
};


/*! A NonlinearSystem together with some TimeDiscretization scheme.
 *
 * This class represents a specific type of first-order ODE, as indicated by
 * \c ODETag.
 *
 * \tparam Matrix the type of matrices occuring in the linearization of the ODE.
 * \tparam Vector the type of the solution vector of the ODE.
 * \tparam ODETag a tag indicating the type of ODE.
 * \tparam NLTag  a tag indicating the method used for resolving nonlinearities.
 */
template<typename Matrix, typename Vector, ODESystemTag ODETag, NonlinearSolverTag NLTag>
class TimeDiscretizedODESystem;


/*! Time discretized first order implicit quasi-linear ODE;
 *  to be solved using the Newton-Raphson method for resolving nonlinearities.
 *
 * \tparam Matrix the type of matrices occuring in the linearization of the ODE.
 * \tparam Vector the type of the solution vector of the ODE.
 *
 * \see ODESystemTag::FirstOrderImplicitQuasilinear
 */
template<typename Matrix, typename Vector>
class TimeDiscretizedODESystem<Matrix, Vector,
                               ODESystemTag::FirstOrderImplicitQuasilinear,
                               NonlinearSolverTag::Newton> final
        : public TimeDiscretizedODESystemBase<Matrix, Vector, NonlinearSolverTag::Newton>
{
public:
    //! A tag indicating the type of ODE.
    static const ODESystemTag ODETag = ODESystemTag::FirstOrderImplicitQuasilinear;

    //! The type of ODE.
    using ODE = ODESystem<Matrix, Vector, ODETag, NonlinearSolverTag::Newton>;
    //! The auxiliary class that computes the matrix/vector used by the nonlinear solver.
    using MatTrans = MatrixTranslator<Matrix, Vector, ODETag>;
    //! A shortcut for a general time discretization scheme
    using TimeDisc = TimeDiscretization<Vector>;

    /*! Constructs a new instance.
     *
     * \param ode the ODE to be wrapped.
     * \param time_discretization the time discretization to be used.
     */
    explicit
    TimeDiscretizedODESystem(ODE& ode, TimeDisc& time_discretization)
        : _ode(ode)
        , _time_disc(time_discretization)
        , _mat_trans(createMatrixTranslator<Matrix, Vector, ODETag>(time_discretization))
        , _Jac(ode.getNumEquations(), ode.getNumEquations())
        , _M  (ode.getNumEquations(), ode.getNumEquations())
        , _K  (ode.getNumEquations(), ode.getNumEquations())
        , _b  (ode.getNumEquations())
    {}

    void assembleResidualNewton(const Vector &x_new_timestep) override
    {
        auto const  t      = _time_disc.getCurrentTime();
        auto const& x_curr = _time_disc.getCurrentX(x_new_timestep);

        _M.setZero();
        _K.setZero();
        _b.setZero();

        _ode.assemble(t, x_curr, _M, _K, _b);
    }

    void assembleJacobian(const Vector &x_new_timestep) override
    {
        namespace BLAS = MathLib::BLAS;

        auto const  t        = _time_disc.getCurrentTime();
        auto const& x_curr   = _time_disc.getCurrentX(x_new_timestep);
        auto const  dxdot_dx = _time_disc.getNewXWeight();
        auto const  dx_dx    = _time_disc.getDxDx();
        _time_disc.getXdot(x_new_timestep, _xdot);

        _Jac.setZero();

        _ode.assembleJacobian(t, x_curr, _xdot,
                              dxdot_dx, _M, dx_dx, _K,
                              _Jac);
    }

    void getResidual(Vector const& x_new_timestep, Vector& res) const override
    {
        // TODO Maybe the duplicate calculation of xdot here and in assembleJacobian
        //      can be optimuized. However, that would make the interface a bit more
        //      fragile.
        _time_disc.getXdot(x_new_timestep, _xdot);

        _mat_trans->computeResidual(_M, _K, _b, x_new_timestep, _xdot, res);
    }

    void getJacobian(Matrix& Jac) const override
    {
        _mat_trans->computeJacobian(_Jac, Jac);
    }

    void applyKnownSolutionsNewton(Matrix& Jac, Vector& res,
                                    Vector& minus_delta_x) override
    {
        (void) Jac; (void) res; (void) minus_delta_x;
        // TODO implement
    }

    bool isLinear() const override
    {
        return _time_disc.isLinearTimeDisc() || _ode.isLinear();
    }

    void pushMatrices() const override
    {
        _mat_trans->pushMatrices(_M, _K, _b);
    }

    TimeDisc& getTimeDiscretization() override {
        return _time_disc;
    }

    MathLib::MatrixSpecifications getMatrixSpecifications() const override {
        return _ode.getMatrixSpecifications();
    }
    std::size_t getNumEquations() const override {
        return _ode.getNumEquations();
    }

private:
    ODE& _ode;            //!< ode the ODE being wrapped
    TimeDisc& _time_disc; //!< the time discretization to being used

    //! the object used to compute the matrix/vector for the nonlinear solver
    std::unique_ptr<MatTrans> _mat_trans;

    Matrix _Jac;  //!< the Jacobian of the residual
    Matrix _M;    //!< Matrix \f$ M \f$.
    Matrix _K;    //!< Matrix \f$ K \f$.
    Vector _b;    //!< Matrix \f$ b \f$.

    // mutable because xdot is computed in the getResidual() method
    mutable Vector _xdot; //!< Used to cache \f$ \dot x \f$. \todo Save some memory.
};


/*! Time discretized first order implicit quasi-linear ODE;
 *  to be solved using the Picard fixpoint iteration method for resolving nonlinearities.
 *
 * \tparam Matrix the type of matrices occuring in the linearization of the ODE.
 * \tparam Vector the type of the solution vector of the ODE.
 *
 * \see ODESystemTag::FirstOrderImplicitQuasilinear
 */
template<typename Matrix, typename Vector>
class TimeDiscretizedODESystem<Matrix, Vector,
                               ODESystemTag::FirstOrderImplicitQuasilinear,
                               NonlinearSolverTag::Picard> final
        : public TimeDiscretizedODESystemBase<Matrix, Vector, NonlinearSolverTag::Picard>
{
public:
    //! A tag indicating the type of ODE.
    static const ODESystemTag ODETag = ODESystemTag::FirstOrderImplicitQuasilinear;

    //! The type of ODE.
    using ODE = ODESystem<Matrix, Vector, ODETag, NonlinearSolverTag::Picard>;
    //! The auxiliary class that computes the matrix/vector used by the nonlinear solver.
    using MatTrans = MatrixTranslator<Matrix, Vector, ODETag>;
    //! A shortcut for a general time discretization scheme
    using TimeDisc = TimeDiscretization<Vector>;

    /*! Constructs a new instance.
     *
     * \param ode the ODE to be wrapped.
     * \param time_discretization the time discretization to be used.
     */
    explicit
    TimeDiscretizedODESystem(ODE& ode, TimeDisc& time_discretization)
        : _ode(ode)
        , _time_disc(time_discretization)
        , _mat_trans(createMatrixTranslator<Matrix, Vector, ODETag>(time_discretization))
        , _M(ode.getNumEquations(), ode.getNumEquations())
        , _K(ode.getNumEquations(), ode.getNumEquations())
        , _b(ode.getNumEquations())
    {}

    void assembleMatricesPicard(const Vector &x_new_timestep) override
    {
        auto const  t      = _time_disc.getCurrentTime();
        auto const& x_curr = _time_disc.getCurrentX(x_new_timestep);

        _M.setZero();
        _K.setZero();
        _b.setZero();

        _ode.assemble(t, x_curr, _M, _K, _b);
    }

    void getA(Matrix& A) const override
    {
        _mat_trans->computeA(_M, _K, A);
    }

    void getRhs(Vector& rhs) const override
    {
        _mat_trans->computeRhs(_M, _K, _b, rhs);
    }

    void applyKnownSolutionsPicard(Matrix& A, Vector& rhs, Vector& x) override
    {
        auto const* known_solutions = _ode.getKnownSolutions();

        if (known_solutions) {
            for (auto const& bc : *known_solutions) {
                // TODO maybe it would be faster to apply all at once
                MathLib::applyKnownSolution(A, rhs, x, bc.global_ids, bc.values);
            }
        }
    }

    bool isLinear() const override
    {
        return _time_disc.isLinearTimeDisc() || _ode.isLinear();
    }

    void pushMatrices() const override
    {
        _mat_trans->pushMatrices(_M, _K, _b);
    }

    TimeDisc& getTimeDiscretization() override {
        return _time_disc;
    }

    MathLib::MatrixSpecifications getMatrixSpecifications() const override {
        return _ode.getMatrixSpecifications();
    }
    std::size_t getNumEquations() const override {
        return _ode.getNumEquations();
    }

private:
    ODE& _ode;            //!< ode the ODE being wrapped
    TimeDisc& _time_disc; //!< the time discretization to being used

    //! the object used to compute the matrix/vector for the nonlinear solver
    std::unique_ptr<MatTrans> _mat_trans;

    Matrix _M;    //!< Matrix \f$ M \f$.
    Matrix _K;    //!< Matrix \f$ K \f$.
    Vector _b;    //!< Matrix \f$ b \f$.
};

//! @}

}

#endif // NUMLIB_TIMEDISCRETIZEDODESYSTEM_H
