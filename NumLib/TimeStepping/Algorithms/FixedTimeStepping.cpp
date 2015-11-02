/**
 * \author Norihiro Watanabe
 * \date   2012-08-03
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "FixedTimeStepping.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <limits>
#include <cassert>

#include "BaseLib/ConfigTreeNew.h"
#include "logog/include/logog.hpp"

namespace NumLib
{

FixedTimeStepping::FixedTimeStepping(double t0, double tn, const std::vector<double> &vec_all_dt)
: _t_initial(t0), _t_end(computeEnd(t0, tn, vec_all_dt)), _dt_vector(vec_all_dt), _ts_prev(t0), _ts_current(t0)
{}

FixedTimeStepping::FixedTimeStepping(double t0, double tn, double dt)
: _t_initial(t0), _t_end(tn), _dt_vector(static_cast<std::size_t>(std::ceil((tn-t0)/dt)), dt), _ts_prev(t0), _ts_current(t0)
{}

std::unique_ptr<FixedTimeStepping>
FixedTimeStepping::newInstance(BaseLib::ConfigTreeNew const& config)
{
    config.checkConfParam("type", "FixedTimeStepping");

    auto const t_initial = config.getConfParam<double>("t_initial");
    auto const t_end     = config.getConfParam<double>("t_end");
    auto const dt        = config.getConfParam<double>("dt");

    std::vector<double> timesteps;
    double te = *t_initial;
    unsigned c = 0;
    double dt = 0.0;

    auto const range = delta_ts->equal_range("pair");
    if (range.first == range.second) {
        ERR("no timesteps have been given");
        return nullptr;
    }
    for (auto it=range.first; it!=range.second; ++it)
    {
        auto count   = it->second.get_optional<std::size_t>("count");
        auto delta_t = it->second.get_optional<double>("delta_t");

        if (count && delta_t) {
            if (*count == 0) {
                ERR("<count> is zero.");
                return nullptr;
            }
            if (*delta_t < 0.0) {
                ERR("timestep <delta_t> is smaller than zero.");
                return nullptr;
            }
            c = *count; dt = *delta_t;
            timesteps.resize(timesteps.size() + c, dt);
        } else {
            ERR("<count> or <delta_t> missing in <pair>");
            return nullptr;
        }

        te += c * dt;

        if (te > *t_end) break;
    }

    // append last timestep until t_end is reached
    if (te < *t_end)
    {
        c = std::ceil((*t_end - te) / dt);
        timesteps.resize(timesteps.size() + c, dt);
        te += c*dt;
    }

    return std::unique_ptr<FixedTimeStepping>(new FixedTimeStepping(*t_initial, *t_end, timesteps));
}

const TimeStep FixedTimeStepping::getTimeStep() const
{
    return _ts_current;
}

bool FixedTimeStepping::next()
{
    // check if last time step
    if (_ts_current.steps() == _dt_vector.size()
        || std::abs(_ts_current.current()-_t_end) < std::numeric_limits<double>::epsilon())
        return false;

    // confirm current time and move to the next if accepted
    if (accepted())
        _ts_prev = _ts_current;

    // prepare the next time step info
    _ts_current = _ts_prev;
    double dt = _dt_vector[_ts_prev.steps()];
    if (_ts_prev.current() + dt > _t_end) // upper bound by t_end
        dt = _t_end - _ts_prev.current();
    _ts_current += dt;

    return true;
}

double FixedTimeStepping::computeEnd(double t_initial, double t_end, const std::vector<double> &dt_vector)
{
    double t_sum = t_initial + std::accumulate(dt_vector.begin(), dt_vector.end(), 0.);
    return std::min(t_end, t_sum);
}

} //NumLib
