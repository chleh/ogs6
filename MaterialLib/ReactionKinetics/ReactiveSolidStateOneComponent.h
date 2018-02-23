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
class ReactiveSolidStateOneComponent : public ReactiveSolidState
{
public:
    explicit ReactiveSolidStateOneComponent(
        std::vector<double> const& solid_density)
        : _solid_density{MathLib::toVector(solid_density)}
    {
        assert(_solid_density.size() == 1);
    }

    ReactiveSolidStateOneComponent& operator+=(
        Eigen::VectorXd const& rhs) override
    {
        _solid_density += rhs;
        return *this;
    }

    ReactiveSolidStateOneComponent& operator+=(
        Eigen::Ref<Eigen::VectorXd> const& rhs) override
    {
        _solid_density += rhs;
        return *this;
    }

    Eigen::VectorXd const& conversion() const override
    {
        return _solid_density;
    }

    Eigen::VectorXd& conversion() override { return _solid_density; }

private:
    Eigen::VectorXd _solid_density;
};

}  // namespace MaterialLib
