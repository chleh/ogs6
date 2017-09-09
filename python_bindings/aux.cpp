#include "reflect-lib/pybind.h"

// empty module used to register auxiliary bindings
REFLECT_LIB_PYTHON_MODULE(AUX, m)
{
    (void)m;
}
