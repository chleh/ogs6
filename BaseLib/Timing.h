/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <chrono>
#include <string>

#include "logog/include/logog.hpp"

namespace BaseLib
{

namespace chr = std::chrono;

enum class TimingType { OneShot, Reuse };

template<enum TimingType>
class Timing;

using TimingOneShot = Timing<TimingType::OneShot>;

template<>
class Timing<TimingType::OneShot>
{
public:
    using Clock = chr::steady_clock;
    using Time = chr::time_point<Clock>;

    Timing(std::string const& label)
        : _start_time{Clock::now()}
        , _label{label}
    {}

    void stop() const
    {
        chr::duration<double> dt = Clock::now() - _start_time;
        INFO("[timing] %s: %f", _label.c_str(), dt.count());
    }

private:
    Time _start_time;
    std::string _label;
};

}
