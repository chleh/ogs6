/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "ReactiveSolidState.h"

namespace MaterialLib
{
class ReactiveSolidStateTwoComponents : public ReactiveSolidState
{
public:
    explicit ReactiveSolidStateTwoComponents(
        std::vector<double> const& conversion)
        : _conversion{MathLib::toVector(conversion)}
    {
        assert(_conversion.size() == 2);
    }

    ReactiveSolidState& operator+=(Eigen::VectorXd const& rhs) override
    {
        _conversion += rhs;
        return *this;
    }
    ReactiveSolidState& operator+=(
        Eigen::Ref<Eigen::VectorXd> const& rhs) override
    {
        _conversion += rhs;
        return *this;
    }

    Eigen::VectorXd const& conversion() const override { return _conversion; }
    Eigen::VectorXd& conversion() override { return _conversion; }

private:
    Eigen::VectorXd _conversion;
};

}  // namespace MaterialLib
