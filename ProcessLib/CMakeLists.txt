include(${PROJECT_SOURCE_DIR}/scripts/cmake/OGSEnabledElements.cmake)

# Collect the process libraries in interface library
add_library(Processes INTERFACE)

if (0)  # TODO [DUNE] re-enable
foreach(process ${ProcessesList})
    if(OGS_BUILD_PROCESS_${process})
        add_subdirectory(${process})
        target_link_libraries(Processes INTERFACE ${process})
    endif()
endforeach()
endif()

add_subdirectory(SmallDeformationDUNE)
target_link_libraries(Processes INTERFACE SmallDeformationDUNE)

add_subdirectory(BoundaryCondition)

APPEND_SOURCE_FILES(SOURCES)
APPEND_SOURCE_FILES(SOURCES SurfaceFlux)
APPEND_SOURCE_FILES(SOURCES Output)
APPEND_SOURCE_FILES(SOURCES Parameter)
APPEND_SOURCE_FILES(SOURCES SourceTerms)
APPEND_SOURCE_FILES(SOURCES Utils)

add_library(ProcessLib ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS ProcessLib LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

target_link_libraries(ProcessLib
    PUBLIC BaseLib MaterialLib MathLib MeshLib NumLib logog MeshGeoToolsLib
    PRIVATE BoundaryCondition
)

if(OGS_USE_PYTHON)
    add_subdirectory(BoundaryCondition/Python)
    target_link_libraries(ProcessLib
        PUBLIC ProcessLibBoundaryConditionPython)
endif()

if(OGS_INSITU)
    target_link_libraries(ProcessLib InSituLib)
endif()

if(OGS_USE_PCH)
    cotire(ProcessLib)
endif()
