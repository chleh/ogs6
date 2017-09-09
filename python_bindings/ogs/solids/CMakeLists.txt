reflect_lib_add_python_module(ogs.solids solids.cpp)
reflect_lib_target_link_libraries(ogs.solids PRIVATE reflection PUBLIC MaterialLib)
add_dependencies(
    REFLECT_LIB_INTERNAL__ogs__solids
    REFLECT_LIB_INTERNAL__ogs__ProcessLib
    REFLECT_LIB_INTERNAL__ogs)
