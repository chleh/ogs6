/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * The code of this file is used to decouple the evaluation of matrix elements from the rest of OGS6,
 * not all of OGS6 has to be recompiled every time a small change is done.
 */

#include "TESFEM-notpl.h"

const double GAS_CONST = 8.3144621;

template <typename T>
T square(const T& v)
{
	return v * v;
}



static inline double fluid_cp(const double p, const double T, const double x)
{
	return 888.888;
}

static inline double fluid_rho(const double p, const double T, const double x)
{
	return 888.888;
}

static inline double solid_cp(const double rho_SR)
{
	return 888.888;
}


namespace ProcessLib
{

namespace TES
{

Eigen::Matrix3d
LADataNoTpl::
getMassCoeffMatrix()
{
	const double rho_GR = fluid_rho(_p, _T, _x);

	const double cpG   = fluid_cp(_p, _T, _x);
	const double cpS   = solid_cp(_rho_SR);

	const double M_pp = _poro/_p * rho_GR;
	const double M_pT = -_poro/_T *  rho_GR;

	double dxn_dxm = _M_inert * _M_react // 0 is inert, 1 is reactive
					 / square(_M_inert * _x + _M_react * (1.0 - _x));

	const double M_px = (_M_react-_M_inert) * _p / (GAS_CONST * _T) * dxn_dxm * _poro;

	const double M_Tp = -_poro;

	const double M_TT = _poro * rho_GR * cpG + (1.0-_poro) * _rho_SR * cpS;
	const double M_Tx = 0.0;

	const double M_xp = 0.0;
	const double M_xT = 0.0;
	const double M_xx = _poro * rho_GR;


	Eigen::Matrix3d M;
	M << M_pp, M_pT, M_px,
		 M_Tp, M_TT, M_Tx,
		 M_xp, M_xT, M_xx;

	return M;
}

} // namespace TES

} // namespace ProcessLib

