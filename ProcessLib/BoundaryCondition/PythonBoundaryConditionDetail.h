/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

namespace ProcessLib
{
struct PyNotOverridden
{
};
class PyBoundaryCondition
{
public:
    virtual std::pair<bool, double> getDirichletBCValue(
        double /*t*/, std::array<double, 3> /*x*/, std::size_t /*node_id*/,
        std::vector<double> const& /*primary_variables*/) const
    {
        _overridden_essential = false;
        return {false, std::numeric_limits<double>::quiet_NaN()};
    }

    virtual std::pair<bool, double> getFlux(
        double /*t*/,
        std::array<double, 3> /*x*/,
        Eigen::VectorXd const& /*primary_variables*/) const
    {
        _overridden_natural = false;
        return {false, std::numeric_limits<double>::quiet_NaN()};
    }

    bool isOverriddenEssential() const { return _overridden_essential; }
    bool isOverriddenNatural() const { return _overridden_natural; }

private:
    mutable bool _overridden_essential = true;
    mutable bool _overridden_natural = true;
};

class PyBoundaryConditionImpl : public PyBoundaryCondition
{
public:
    using PyBoundaryCondition::PyBoundaryCondition;

    std::pair<bool, double> getDirichletBCValue(
        double t, std::array<double, 3> x, std::size_t node_id,
        std::vector<double> const& primary_variables) const override
    {
        using Ret = std::pair<bool, double>;
        PYBIND11_OVERLOAD(Ret, PyBoundaryCondition, getDirichletBCValue, t, x,
                          node_id, primary_variables);
    }

    virtual std::pair<bool, double> getFlux(
        double t, std::array<double, 3> x,
        Eigen::VectorXd const& primary_variables) const
    {
        using Ret = std::pair<bool, double>;
        PYBIND11_OVERLOAD(Ret, PyBoundaryCondition, getFlux, t, x,
                          primary_variables);
    }
};
}  // namespace ProcessLib
