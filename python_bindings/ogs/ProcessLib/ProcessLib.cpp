#include "reflect-lib/pybind.h"

#include "ProcessLib/Parameter/ConstantParameter.h"

REFLECT_LIB_PYTHON_MODULE(ogs__ProcessLib, module)
{
    using namespace ProcessLib;
    namespace py = pybind11;

    reflect_lib::Module m(module);

    m.bind<ParameterBase>();
    m.bind<Parameter<double>>();
    m.bind<ConstantParameter<double>>()
        .def(py::init<std::string const&, double>())
        .def(py::init<std::string const&, std::vector<double>>());
    m.bind<SpatialPosition>();
}
