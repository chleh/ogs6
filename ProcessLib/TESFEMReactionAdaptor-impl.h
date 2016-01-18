/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <cassert>

#include <logog/include/logog.hpp>

#include "TESFEMReactionAdaptor.h"
#include "MathLib/Nonlinear/Root1D.h"

#include "MaterialsLib/adsorption/adsorption.h"
#include "MaterialsLib/adsorption/reaction_inert.h"
#include "MaterialsLib/adsorption/reaction_sinusoidal.h"

namespace ProcessLib
{
namespace TES
{

template<typename Traits>
std::unique_ptr<TESFEMReactionAdaptor<Traits> >
TESFEMReactionAdaptor<Traits>::
newInstance(LADataNoTpl<Traits>& data)
{
    auto const* ads = data._AP->_reaction_system.get();
    if (dynamic_cast<Ads::Adsorption const*>(ads) != nullptr) {
        return std::unique_ptr<TESFEMReactionAdaptor<Traits> >(
                    new TESFEMReactionAdaptorAdsorption<Traits>(data)
                    );
    } else if (dynamic_cast<Ads::ReactionInert const*>(ads) != nullptr) {
        return std::unique_ptr<TESFEMReactionAdaptor<Traits> >(
                    new TESFEMReactionAdaptorInert<Traits>(data)
                    );
    } else if (dynamic_cast<Ads::ReactionSinusoidal const*>(ads) != nullptr) {
        return std::unique_ptr<TESFEMReactionAdaptor<Traits> >(
                    new TESFEMReactionAdaptorSinusoidal<Traits>(data)
                    );
    } else if (dynamic_cast<Ads::ReactionCaOH2 const*>(ads) != nullptr) {
        return std::unique_ptr<TESFEMReactionAdaptor<Traits> >(
                    new TESFEMReactionAdaptorCaOH2<Traits>(data)
                    );
    }

    ERR("No suitable TESFEMReactionAdaptor found. Aborting.");
    std::abort();
    return std::unique_ptr<TESFEMReactionAdaptor<Traits> >(nullptr);
}


template<typename Traits>
TESFEMReactionAdaptorAdsorption<Traits>::
TESFEMReactionAdaptorAdsorption(LADataNoTpl<Traits> &data)
    // caution fragile: this relies in this constructor b eing called __after__
    // data._solid_density has been properly set up!
    : _bounds_violation(data._solid_density.size(), false)
    , _data(data)
{
    assert(dynamic_cast<Ads::Adsorption const*>(data._AP->_reaction_system.get()) != nullptr
           && "Reactive system has wrong type.");
    assert(_bounds_violation.size() != 0);
}

template<typename Traits>
void
TESFEMReactionAdaptorAdsorption<Traits>::
initReaction_slowDownUndershootStrategy(const unsigned int_pt)
{
    auto const& AP = _data._AP;
    assert(AP->_number_of_try_of_iteration < 20);

    const double loading = Ads::Adsorption::get_loading(_data._solid_density_prev_ts[int_pt],
                                                        AP->_rho_SR_dry);

    double react_rate_R = AP->_reaction_system->get_reaction_rate(_data._p_V, _data._T, AP->_M_react, loading)
                          * AP->_rho_SR_dry;

    // set reaction rate based on current damping factor
    react_rate_R = (_reaction_damping_factor > 1e-3)
                   ? _reaction_damping_factor * react_rate_R
                   : 0.0;

    if (_data._p_V < 0.01 * Ads::Adsorption::get_equilibrium_vapour_pressure(_data._T)
        && react_rate_R > 0.0)
    {
        react_rate_R = 0.0;
    }
    else if (_data._p_V < 100.0
             || _data._p_V < 0.05 * Ads::Adsorption::get_equilibrium_vapour_pressure(_data._T))
    {
        // use equilibrium reaction for dry regime

        // in the case of zeroth try in zeroth iteration: _p_V and loading are the values
        // at the end of the previous timestep

        const double pV_eq = estimateAdsorptionEquilibrium(_data._p_V, loading);
        // TODO [CL]: it would be more correct to subtract pV from the previous timestep here
        const double delta_pV = pV_eq - _data._p_V;
        const double delta_rhoV = delta_pV * AP->_M_react / Ads::GAS_CONST / _data._T * AP->_poro;
        const double delta_rhoSR = delta_rhoV / (AP->_poro - 1.0);
        double react_rate_R2 = delta_rhoSR/AP->_delta_t;

        if (_bounds_violation[int_pt])
        {
            react_rate_R2 *= 0.5;
        }

        // 0th try: make sure reaction is not slower than allowed by local estimation
        // nth try: make sure reaction is not made faster by local estimation
        if (
            (AP->_number_of_try_of_iteration == 0
             && std::abs(react_rate_R2) > std::abs(react_rate_R))
            ||
            (AP->_number_of_try_of_iteration != 0
             && std::abs(react_rate_R2) < std::abs(react_rate_R))
            )
        {
            react_rate_R = react_rate_R2;
        }
    }

    // smooth out readjustment of reaction rate
    if (AP->_iteration_in_current_timestep > 3)
    {
        if (AP->_iteration_in_current_timestep <= 8)
        {
            // update reaction rate for for five iterations
            const auto N = AP->_iteration_in_current_timestep - 3;

            // take average s.t. does not oscillate so much
            react_rate_R = 1.0 / (1.0+N) * (N*_data._reaction_rate[int_pt] + 1.0 * react_rate_R);
        }
        else
        {
            // afterwards no update anymore
            react_rate_R = _data._reaction_rate[int_pt];
        }
    }

    if (AP->_number_of_try_of_iteration > 0)
    {
        // assert that within tries reaction does not get faster
        // (e.g. due to switch equilibrium reaction <--> kinetic reaction)

        // factor of 0.9*N: in fact, even slow down reaction over tries
        const double r = std::pow(0.9, AP->_number_of_try_of_iteration)
                         *_data._reaction_rate[int_pt];
        if (std::abs(react_rate_R) > std::abs(r)) {
            react_rate_R = r;
        }
    }

    _data._reaction_rate[int_pt] = react_rate_R;
    _data._solid_density[int_pt] = _data._solid_density_prev_ts[int_pt] + react_rate_R * AP->_delta_t;

    _data._qR = _data._reaction_rate[int_pt];
}

template<typename Traits>
double
TESFEMReactionAdaptorAdsorption<Traits>::
estimateAdsorptionEquilibrium(const double p_V0, const double C0) const
{
    auto const &data = this->_data;

    auto f = [&data, p_V0, C0](double pV) -> double
    {
        auto const& AP = data._AP;
        // pV0 := _p_V
        const double C_eq = AP->_reaction_system->get_equilibrium_loading(pV, data._T, AP->_M_react);
        return (pV - p_V0) * AP->_M_react / Ads::GAS_CONST / data._T * AP->_poro
                + (1.0-AP->_poro) * (C_eq - C0) * AP->_rho_SR_dry;
    };

    // range where to search for roots of f
    const double C_eq0 = data._AP->_reaction_system->get_equilibrium_loading(p_V0, data._T, data._AP->_M_react);
    const double limit = (C_eq0 > C0)
                         ? 1e-8
                         : Ads::Adsorption::get_equilibrium_vapour_pressure(data._T);

    // search for roots
    auto rf = MathLib::Nonlinear::makeRegulaFalsi<MathLib::Nonlinear::Pegasus>(f, p_V0, limit);
    rf.step(3);

    // set vapour pressure
    return rf.get_result();
}

template<typename Traits>
bool
TESFEMReactionAdaptorAdsorption<Traits>::
checkBounds(std::vector<double> const& localX,
            std::vector<double> const& localX_pts)
{
    double alpha = 1.0;

    const double min_xmV = 1e-6;
    const std::size_t nnodes = localX.size() / NODAL_DOF;
    const std::size_t xmV_offset = (NODAL_DOF - 1)*nnodes;

    for (std::size_t i=0; i<nnodes; ++i)
    {
        auto const xnew = localX[xmV_offset+i];
        if (xnew < min_xmV)
        {
            auto const xold = localX_pts[i+xmV_offset];
            const auto a = xold / (xold - xnew);
            alpha = std::min(alpha, a);
            _bounds_violation[i] = true;
        }
        else if (xnew > 1.0)
        {
            auto const xold = localX_pts[i+xmV_offset];
            const auto a = xold / (xnew - xold);
            alpha = std::min(alpha, a);
            _bounds_violation[i] = true;
        }
        else
        {
            _bounds_violation[i] = false;
        }
    }

    assert (alpha > 0.0);

    if (alpha != 1.0)
    {
        if (_data._AP->_number_of_try_of_iteration <=2) {
            _reaction_damping_factor *= sqrt(std::min(alpha, 0.5));
        } else {
            _reaction_damping_factor *= std::min(alpha, 0.5);
        }
    }

    return alpha == 1.0;
}

template<typename Traits>
void
TESFEMReactionAdaptorAdsorption<Traits>::
preZerothTryAssemble()
{
    _reaction_damping_factor = std::min(
        std::sqrt(_reaction_damping_factor),
        10.0*_reaction_damping_factor);
}


template<typename Traits>
TESFEMReactionAdaptorInert<Traits>::
TESFEMReactionAdaptorInert(LADataNoTpl<Traits>& data)
    : _data(data)
{
}


template<typename Traits>
TESFEMReactionAdaptorSinusoidal<Traits>::
TESFEMReactionAdaptorSinusoidal(LADataNoTpl<Traits> &data)
    : _data{data}
{
    assert(dynamic_cast<Ads::ReactionSinusoidal const*>(data._AP->_reaction_system.get()) != nullptr
           && "Reactive system has wrong type.");
}


template<typename Traits>
void
TESFEMReactionAdaptorSinusoidal<Traits>::
initReaction(const unsigned int int_pt)
{
    const double t = _data._AP->_current_time;

    // Cf. OGS5
    const double rhoSR0 = 1.0;
    const double rhoTil = 0.1;
    const double omega  = 2.0 * 3.1416;
    const double poro   = _data._AP->_poro;

    _data._solid_density[int_pt] = rhoSR0 + rhoTil * std::sin(omega*t) / (1.0 - poro);
    _data._reaction_rate[int_pt] = rhoTil * omega * cos(omega*t) / (1.0 - poro);
    _data._qR = _data._reaction_rate[int_pt];
}


template<typename Traits>
TESFEMReactionAdaptorCaOH2<Traits>::
TESFEMReactionAdaptorCaOH2(LADataNoTpl<Traits> &data)
    : _data{data}
    , _react{dynamic_cast<Ads::ReactionCaOH2&>(*data._AP->_reaction_system.get())}
{
    _ode_solver = std::move(MathLib::createOdeSolver<1, React>(_react.getOdeSolverConfig()));

    _ode_solver->init();
    _ode_solver->setTolerance(1e-10, 1e-10);

    _ode_solver->setFunction(odeRhs, nullptr, &_react); // TODO: change signature to reference
}

template<typename Traits>
void
TESFEMReactionAdaptorCaOH2<Traits>::
initReaction(const unsigned int int_pt)
{
    if (_data._AP->_iteration_in_current_timestep > 0
        || _data._AP->_number_of_try_of_iteration > 0)
    {
        _data._qR = _data._reaction_rate[int_pt];
        return;
    }


    // TODO: double check!
	// const double xv_NR  = SolidProp->non_reactive_solid_volume_fraction;
	// const double rho_NR = SolidProp->non_reactive_solid_density;
	const double xv_NR = 0.0;
	const double rho_NR = 0.0;


	const double t0 = 0.0;
	const double y0 = (_data._solid_density_prev_ts[int_pt] - xv_NR * rho_NR) / (1.0-xv_NR);

	const double t_end = _data._AP->_delta_t;

	_react.update_param(_data._T, _data._p, _data._vapour_mass_fraction,
						_data._solid_density_prev_ts[int_pt]);

	_ode_solver->setIC(t0, { y0 });
	_ode_solver->preSolve();
	_ode_solver->solve(t_end);

	const double time_reached = _ode_solver->getTime();
	assert(std::abs(t_end - time_reached) < std::numeric_limits<double>::epsilon());

	const double y_new     = _ode_solver->getSolution()[0];
	const double y_dot_new = _ode_solver->getYDot(t_end, { y_new })[0];

	double rho_react;

	//cut off when limits are reached
	if ( y_new < _react.rho_low )
		rho_react = _react.rho_low;
	else if ( y_new > _react.rho_up ) //{
		rho_react = _react.rho_up;
	else
		rho_react = y_new;

    _data._solid_density[int_pt] = (1.0-xv_NR) * rho_react + xv_NR * rho_NR;
    _data._reaction_rate[int_pt] = y_dot_new * (1.0-xv_NR);
    _data._qR = _data._reaction_rate[int_pt];
}


// TODO: there could be a general implementation for this in the OdeSolver class
// general implementation should take an object and a member function pointer
// with suitable signature.
template<typename Traits>
bool
TESFEMReactionAdaptorCaOH2<Traits>::
odeRhs(const double t,
       BaseLib::ArrayRef<const double, 1> const y,
       BaseLib::ArrayRef<double, 1> ydot,
       Ads::ReactionCaOH2& reaction)
{
    reaction.eval(t, y, ydot);
    return true;
}


}
}
