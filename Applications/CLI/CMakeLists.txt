
add_executable(ogs ogs.cpp)

target_link_libraries(ogs
    PRIVATE BaseLib ApplicationsLib
)

if(OGS_USE_PYTHON)
    message(STATUS "OGS uses Python!!! ${PYTHON_LIBRARIES}")
    target_link_libraries(ogs PRIVATE pybind11::embed)
    target_link_libraries(ogs PRIVATE ${PYTHON_LIBRARIES})
    # target_include_directories(ogs PRIVATE ${PYTHON_INCLUDE_DIRS})

    add_library(ogs_python_bindings STATIC ogs_python_bindings.cpp)
    target_link_libraries(ogs_python_bindings PRIVATE pybind11::embed)
    target_include_directories(ogs_python_bindings PRIVATE ${PYTHON_INCLUDE_DIRS})

    target_link_libraries(ogs PRIVATE ogs_python_bindings)

    if(BUILD_SHARED_LIBS)
        target_compile_definitions(ogs_python_bindings OGS_BUILD_SHARED_LIBS)
    endif()
endif()

if(OGS_USE_PETSC)
    target_link_libraries(ogs PRIVATE ${PETSC_LIBRARIES})
endif()

if(OGS_INSITU)
    target_link_libraries(ogs PRIVATE InSituLib)
endif()

if(OGS_USE_PCH)
    cotire(ogs)
endif()

####################
### Tests ##########
####################
add_test(NAME ogs_no_args COMMAND ogs)
set_tests_properties(ogs_no_args PROPERTIES WILL_FAIL TRUE)

####################
### Installation ###
####################
install(TARGETS ogs RUNTIME DESTINATION bin COMPONENT ogs_cli)

set(CPACK_PACKAGE_EXECUTABLES ${CPACK_PACKAGE_EXECUTABLES} "ogs" "OGS Simulator")
cpack_add_component(ogs_cli
    DISPLAY_NAME "OGS THMC Simulator"
    DESCRIPTION "The command line interface for OpenGeoSys."
    GROUP Applications
)
