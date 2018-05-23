#pragma once

#ifdef OGS_USE_PYTHON

#include <pybind11/embed.h>

namespace ApplicationsLib
{
pybind11::scoped_interpreter setupEmbeddedPython();

}  // namespace ApplicationsLib

#endif  // OGS_USE_PYTHON
