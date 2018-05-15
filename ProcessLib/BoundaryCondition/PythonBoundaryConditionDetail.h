/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <pybind11/pybind11.h>

namespace
{
struct PyNotOverridden
{
};
class PyBoundaryCondition
{
public:
    virtual std::pair<bool, double> getDirichletBCValue(
        double /*t*/, std::array<double, 3> /*x*/) const
    {
        throw PyNotOverridden{};
    }

    virtual std::pair<bool, double> getFlux(double /*t*/) const
    {
        throw PyNotOverridden{};
    }
};

class PyBoundaryConditionImpl : public PyBoundaryCondition
{
public:
    using PyBoundaryCondition::PyBoundaryCondition;

    std::pair<bool, double> getDirichletBCValue(
        double t, std::array<double, 3> x) const override
    {
        using Ret = std::pair<bool, double>;
        PYBIND11_OVERLOAD(Ret, PyBoundaryCondition, getDirichletBCValue, t, x);
    }

    virtual std::pair<bool, double> getFlux(double t) const
    {
        using Ret = std::pair<bool, double>;
        PYBIND11_OVERLOAD(Ret, PyBoundaryCondition, getFlux, t);
    }
};
}  // namespace
