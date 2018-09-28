/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "PythonSourceTermModule.h"

#include <pybind11/stl.h>

#include "PythonSourceTermPythonSideInterface.h"

namespace ProcessLib
{
//! Trampoline class allowing methods of class
//! PythonSourceTermPythonSideInterface to be overridden on the Python
//! side. Cf. https://pybind11.readthedocs.io/en/stable/advanced/classes.html
class PythonSourceTermPythonSideInterfaceTrampoline
    : public PythonSourceTermPythonSideInterface
{
public:
    using PythonSourceTermPythonSideInterface::
        PythonSourceTermPythonSideInterface;

    std::pair<double, std::array<double, 3>> getFlux(
        double t, std::array<double, 3> x,
        std::vector<double> const& primary_variables) const override
    {
        using Ret = std::pair<double, std::array<double, 3>>;
        PYBIND11_OVERLOAD(Ret, PythonSourceTermPythonSideInterface,
                          getFlux, t, x, primary_variables);
    }
};

void pythonBindSourceTerm(pybind11::module& m_st)
{
    namespace py = pybind11;

    py::class_<PythonSourceTermPythonSideInterface,
               PythonSourceTermPythonSideInterfaceTrampoline>
        pybc(m_st, "SourceTerm");

    pybc.def(py::init());

    pybc.def("getFlux", &PythonSourceTermPythonSideInterface::getFlux);
}

}  // namespace ProcessLib
