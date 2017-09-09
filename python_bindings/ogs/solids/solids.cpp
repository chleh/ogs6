#include "reflect-lib/pybind.h"

#include "MaterialLib/SolidModels/LinearElasticIsotropic.h"

REFLECT_LIB_PYTHON_MODULE(ogs__solids, module)
{
    using namespace MaterialLib::Solids;
    namespace py = pybind11;

    reflect_lib::Module m(module);

    m.bind<MechanicsBase<2>>();
    m.bind<MechanicsBase<3>>();

    m.bind<LinearElasticIsotropic<2>>();
    m.bind<LinearElasticIsotropic<3>>();

    // Using that manual "trick" MaterialProperties cannot be pickled!
    m.bind<LinearElasticIsotropic<2>::MaterialProperties>().def(
        py::init<ProcessLib::Parameter<double> const&,
                 ProcessLib::Parameter<double> const&>());
    m.bind<LinearElasticIsotropic<3>::MaterialProperties>().def(
        py::init<ProcessLib::Parameter<double> const&,
                 ProcessLib::Parameter<double> const&>());
}
