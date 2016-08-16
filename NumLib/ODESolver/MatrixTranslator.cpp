/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MatrixTranslator.h"

#include "MathLib/LinAlg/LinAlg.h"

namespace NumLib
{
void MatrixTranslatorGeneral<ODESystemTag::FirstOrderImplicitQuasilinear>::
    computeA(GlobalMatrix const& M, GlobalMatrix const& K,
             GlobalMatrix& A) const
{
    namespace LinAlg = MathLib::LinAlg;

    auto const dxdot_dx = _time_disc.getNewXWeight();

    // A = M * dxdot_dx + K
    LinAlg::copy(M, A);
    LinAlg::aypx(A, dxdot_dx, K);
}

void MatrixTranslatorGeneral<ODESystemTag::FirstOrderImplicitQuasilinear>::
    computeRhs(const GlobalMatrix& M, const GlobalMatrix& /*K*/,
               const GlobalVector& b, GlobalVector& rhs) const
{
    namespace LinAlg = MathLib::LinAlg;

    auto& tmp = NumLib::GlobalVectorProvider::provider.getVector(_tmp_id);
    _time_disc.getWeightedOldX(tmp);

    // rhs = M * weighted_old_x + b
    LinAlg::matMultAdd(M, tmp, b, rhs);

    NumLib::GlobalVectorProvider::provider.releaseVector(tmp);
}

void MatrixTranslatorGeneral<ODESystemTag::FirstOrderImplicitQuasilinear>::
    computeResidual(GlobalMatrix const& M, GlobalMatrix const& K,
                    GlobalVector const& b, GlobalVector const& x_new_timestep,
                    GlobalVector const& xdot, GlobalVector& res) const
{
    namespace LinAlg = MathLib::LinAlg;

    auto const& x_curr = _time_disc.getCurrentX(x_new_timestep);

    // res = M * x_dot + K * x_curr - b
    LinAlg::matMult(M, xdot, res);  // the local vector x_dot seems to be
                                  // necessary because of this multiplication
    LinAlg::matMultAdd(K, x_curr, res, res);
    LinAlg::axpy(res, -1.0, b);
}

void MatrixTranslatorGeneral<ODESystemTag::FirstOrderImplicitQuasilinear>::
    computeJacobian(GlobalMatrix const& Jac_in, GlobalMatrix& Jac_out) const
{
    namespace LinAlg = MathLib::LinAlg;

    LinAlg::copy(Jac_in, Jac_out);
}

void MatrixTranslatorForwardEuler<ODESystemTag::FirstOrderImplicitQuasilinear>::
    computeA(GlobalMatrix const& M, GlobalMatrix const& /*K*/,
             GlobalMatrix& A) const
{
    namespace LinAlg = MathLib::LinAlg;

    auto const dxdot_dx = _fwd_euler.getNewXWeight();

    // A = M * dxdot_dx
    LinAlg::copy(M, A);
    LinAlg::scale(A, dxdot_dx);
}

void MatrixTranslatorForwardEuler<ODESystemTag::FirstOrderImplicitQuasilinear>::
    computeRhs(const GlobalMatrix& M, const GlobalMatrix& K,
               const GlobalVector& b, GlobalVector& rhs) const
{
    namespace LinAlg = MathLib::LinAlg;

    auto& tmp = NumLib::GlobalVectorProvider::provider.getVector(_tmp_id);
    _fwd_euler.getWeightedOldX(tmp);

    auto const& x_old = _fwd_euler.getXOld();

    // rhs = b + M * weighted_old_x - K * x_old
    LinAlg::matMult(K, x_old, rhs);        // rhs = K * x_old
    LinAlg::aypx(rhs, -1.0, b);            // rhs = b - K * x_old
    LinAlg::matMultAdd(M, tmp, rhs, rhs);  // rhs += M * weighted_old_x

    NumLib::GlobalVectorProvider::provider.releaseVector(tmp);
}

void MatrixTranslatorForwardEuler<ODESystemTag::FirstOrderImplicitQuasilinear>::
    computeResidual(GlobalMatrix const& M, GlobalMatrix const& K,
                    GlobalVector const& b, GlobalVector const& x_new_timestep,
                    GlobalVector const& xdot, GlobalVector& res) const
{
    namespace LinAlg = MathLib::LinAlg;

    auto const& x_curr = _fwd_euler.getCurrentX(x_new_timestep);

    // res = M * x_dot + K * x_curr - b
    LinAlg::matMult(M, xdot, res);
    LinAlg::matMultAdd(K, x_curr, res, res);
    LinAlg::axpy(res, -1.0, b);
}

void MatrixTranslatorForwardEuler<ODESystemTag::FirstOrderImplicitQuasilinear>::
    computeJacobian(GlobalMatrix const& Jac_in, GlobalMatrix& Jac_out) const
{
    namespace LinAlg = MathLib::LinAlg;

    LinAlg::copy(Jac_in, Jac_out);
}

}  // NumLib
