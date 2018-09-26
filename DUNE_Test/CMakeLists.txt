add_executable("dune-hello" dune-hello.cc)
target_link_dune_default_libraries("dune-hello")

add_test(NAME "dune-hello-test"
    COMMAND "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/dune-hello"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")


add_executable("dune-poisson-eq" dune-poisson-eq.cc)
target_link_dune_default_libraries("dune-poisson-eq")

add_test(NAME "dune-poisson-eq"
    COMMAND
    # TODO [DUNE] search for shell or remove /bin/sh
    "/bin/sh" "-c" "'$<TARGET_FILE:dune-poisson-eq>' '${CMAKE_CURRENT_SOURCE_DIR}/l-shape.msh' && '$<TARGET_FILE:vtkdiff>' --rel 0 --abs 1e-16 -a solution -b solution '${CMAKE_CURRENT_SOURCE_DIR}/poissonequation_result_ref.vtu' '${CMAKE_CURRENT_BINARY_DIR}/poissonequation_result.vtu'"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")


add_executable("dune-tri6-io" dune-tri6-io.cc)
target_link_dune_default_libraries("dune-tri6-io")

add_test(NAME "dune-tri6-io"
    COMMAND
    # TODO [DUNE] search for shell or remove /bin/sh
    "/bin/sh" "-c" "'$<TARGET_FILE:dune-tri6-io>' '${CMAKE_CURRENT_SOURCE_DIR}/square-1x1-tri6.msh' '${CMAKE_CURRENT_BINARY_DIR}/square-1x1-tri6'"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
