/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#if OGS_USE_PYTHON

#include "PythonBoundaryConditionModule.h"

#include "PythonBoundaryConditionDetail.h"

namespace ProcessLib
{
namespace py = pybind11;
void pythonBindBoundaryCondition(py::module& m)
{
    // `m` is a `py::module` which is used to bind functions and classes
    py::class_<::PyBoundaryCondition, ::PyBoundaryConditionImpl> pybc(
        m, "BoundaryCondition");
    pybc.def(py::init());
    pybc.def("getDirichletBCValue",
             &::PyBoundaryCondition::getDirichletBCValue);
    pybc.def("getFlux", &::PyBoundaryCondition::getFlux);
}

}  // namespace ProcessLib

#endif
