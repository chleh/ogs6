APPEND_SOURCE_FILES(SOURCES)

add_library(SmallDeformationDUNE ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS SmallDeformationDUNE LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
target_link_libraries(SmallDeformationDUNE PUBLIC ProcessLib)

# TODO [DUNE] fix: include(Tests.cmake)
