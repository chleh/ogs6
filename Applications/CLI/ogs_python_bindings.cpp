#ifdef OGS_USE_PYTHON

#include <pybind11/embed.h>

#include "ogs_python_bindings.h"

#ifndef OGS_BUILD_SHARED_LIBS
extern "C" PyObject* pybind11_init_impl_OpenGeoSys();

template <typename T>
void mark_used(T p)
{
    volatile T vp = (volatile T)p;
    vp = vp;
}

#endif  // OGS_BUILD_SHARED_LIBS

namespace ApplicationsLib
{
pybind11::scoped_interpreter setupEmbeddedPython()
{
#ifndef OGS_BUILD_SHARED_LIBS
    mark_used(&pybind11_init_impl_OpenGeoSys);
#endif

    return pybind11::scoped_interpreter{};
}

}  // namespace ApplicationsLib

#endif  // OGS_USE_PYTHON
