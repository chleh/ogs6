unset(CMAKE_LIBRARY_OUTPUT_DIRECTORY)
add_compile_options(-ftemplate-backtrace-limit=0)

include_directories("/home/lehmannc/prog/cereal-pybind/src")
add_subdirectory(pybind11)

set(REFLECT_LIB_SOURCE_ROOT "/home/lehmannc/prog/cereal-pybind/src")
include("${REFLECT_LIB_SOURCE_ROOT}/cmake/python_bindings.cmake")


file(GLOB sources LIST_DIRECTORIES false "${REFLECT_LIB_SOURCE_ROOT}/reflect-lib/*.cpp" "${REFLECT_LIB_SOURCE_ROOT}/reflect-lib/*.h")
add_library(reflection ${sources})
target_compile_options(reflection PRIVATE -fPIC)

reflect_lib_add_aux_python_module(AUX aux.cpp)
reflect_lib_target_link_libraries(AUX PRIVATE reflection)

add_subdirectory(ogs)
