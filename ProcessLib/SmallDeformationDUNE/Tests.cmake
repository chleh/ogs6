##### DUNE #####

# Simple test if small deformation works with DUNE
# Stretching in y-direction with a uniform Dirichlet BC
AddTest(
    NAME DUNE_SmallDefLin_square_coarse
    PATH Mechanics/LinearDUNE
    EXECUTABLE_ARGS square_coarse.prj
    REQUIREMENTS NOT OGS_USE_MPI
    TESTER vtkdiff
    ABSTOL 1e-15 RELTOL 0
    DIFF_DATA
    square_coarse_solution.vtu square_coarse_pcs_0_ts_4_t_1.000000.vtu displacement[0] displacement[0]
    square_coarse_solution.vtu square_coarse_pcs_0_ts_4_t_1.000000.vtu displacement[1] displacement[1]
    #
    square_coarse_solution.vtu  square_coarse_pcs_0_ts_4_t_1.000000.vtu sigma[0] sigma[0]
    square_coarse_solution.vtu  square_coarse_pcs_0_ts_4_t_1.000000.vtu sigma[1] sigma[1]
    square_coarse_solution.vtu  square_coarse_pcs_0_ts_4_t_1.000000.vtu sigma[2] sigma[2]
    square_coarse_solution.vtu  square_coarse_pcs_0_ts_4_t_1.000000.vtu sigma[3] sigma[3]
    #
    square_coarse_solution.vtu  square_coarse_pcs_0_ts_4_t_1.000000.vtu epsilon[0] epsilon[0]
    square_coarse_solution.vtu  square_coarse_pcs_0_ts_4_t_1.000000.vtu epsilon[1] epsilon[1]
    square_coarse_solution.vtu  square_coarse_pcs_0_ts_4_t_1.000000.vtu epsilon[2] epsilon[2]
    square_coarse_solution.vtu  square_coarse_pcs_0_ts_4_t_1.000000.vtu epsilon[3] epsilon[3]
    )

# Simple test if global refinement works with DUNE
# Stretching in y-direction with a uniform Dirichlet BC after global mesh refinement
AddTest(
    NAME DUNE_SmallDefLin_square_globally_refined
    PATH Mechanics/LinearDUNE
    EXECUTABLE_ARGS square_globally_refined.prj
    REQUIREMENTS NOT OGS_USE_MPI
    TESTER vtkdiff
    ABSTOL 1e-15 RELTOL 0
    DIFF_DATA
    square_globally_refined_solution.vtu square_globally_refined_pcs_0_ts_4_t_1.000000.vtu displacement[0] displacement[0]
    square_globally_refined_solution.vtu square_globally_refined_pcs_0_ts_4_t_1.000000.vtu displacement[1] displacement[1]
    #
    square_globally_refined_solution.vtu square_globally_refined_pcs_0_ts_4_t_1.000000.vtu sigma[0] sigma[0]
    square_globally_refined_solution.vtu square_globally_refined_pcs_0_ts_4_t_1.000000.vtu sigma[1] sigma[1]
    square_globally_refined_solution.vtu square_globally_refined_pcs_0_ts_4_t_1.000000.vtu sigma[2] sigma[2]
    square_globally_refined_solution.vtu square_globally_refined_pcs_0_ts_4_t_1.000000.vtu sigma[3] sigma[3]
    #
    square_globally_refined_solution.vtu square_globally_refined_pcs_0_ts_4_t_1.000000.vtu epsilon[0] epsilon[0]
    square_globally_refined_solution.vtu square_globally_refined_pcs_0_ts_4_t_1.000000.vtu epsilon[1] epsilon[1]
    square_globally_refined_solution.vtu square_globally_refined_pcs_0_ts_4_t_1.000000.vtu epsilon[2] epsilon[2]
    square_globally_refined_solution.vtu square_globally_refined_pcs_0_ts_4_t_1.000000.vtu epsilon[3] epsilon[3]
    )

# Simple test if Dirichlet values are interpolated correctly during global refinement
# Stretching in y-direction with a non-uniform Dirichlet BC after global mesh refinement
# Hint for verification: It is most important that the y-displacement at the top boundary
# matches the prescribed linear slope exactly!
AddTest(
    NAME DUNE_SmallDefLin_square_globally_refined_gradient
    PATH Mechanics/LinearDUNE
    EXECUTABLE_ARGS square_globally_refined_gradient.prj
    REQUIREMENTS NOT OGS_USE_MPI
    TESTER vtkdiff
    ABSTOL 1e-15 RELTOL 0
    DIFF_DATA
    square_globally_refined_gradient_solution.vtu square_globally_refined_gradient_pcs_0_ts_4_t_1.000000.vtu displacement[0] displacement[0]
    square_globally_refined_gradient_solution.vtu square_globally_refined_gradient_pcs_0_ts_4_t_1.000000.vtu displacement[1] displacement[1]
    #
    square_globally_refined_gradient_solution.vtu square_globally_refined_gradient_pcs_0_ts_4_t_1.000000.vtu sigma[0] sigma[0]
    square_globally_refined_gradient_solution.vtu square_globally_refined_gradient_pcs_0_ts_4_t_1.000000.vtu sigma[1] sigma[1]
    square_globally_refined_gradient_solution.vtu square_globally_refined_gradient_pcs_0_ts_4_t_1.000000.vtu sigma[2] sigma[2]
    square_globally_refined_gradient_solution.vtu square_globally_refined_gradient_pcs_0_ts_4_t_1.000000.vtu sigma[3] sigma[3]
    #
    square_globally_refined_gradient_solution.vtu square_globally_refined_gradient_pcs_0_ts_4_t_1.000000.vtu epsilon[0] epsilon[0]
    square_globally_refined_gradient_solution.vtu square_globally_refined_gradient_pcs_0_ts_4_t_1.000000.vtu epsilon[1] epsilon[1]
    square_globally_refined_gradient_solution.vtu square_globally_refined_gradient_pcs_0_ts_4_t_1.000000.vtu epsilon[2] epsilon[2]
    square_globally_refined_gradient_solution.vtu square_globally_refined_gradient_pcs_0_ts_4_t_1.000000.vtu epsilon[3] epsilon[3]
    )

# Simple test for the error estimator
# Stretching in y-direction with a non-uniform Dirichlet BC and local mesh refinement
# This test also verifies that the mesh refinement procedure comes to an end and does
# not run forever.
# Hint for verification: It is most important that the y-displacement at the top boundary
# matches the prescribed linear slope exactly!
AddTest(
    NAME DUNE_SmallDefLin_square_error_estimator
    PATH Mechanics/LinearDUNE
    EXECUTABLE_ARGS square_error_estimator.prj
    REQUIREMENTS NOT OGS_USE_MPI
    TESTER vtkdiff
    ABSTOL 1e-15 RELTOL 0
    DIFF_DATA
    square_error_estimator_solution.vtu square_error_estimator_pcs_0_ts_4_t_1.000000.vtu displacement[0] displacement[0]
    square_error_estimator_solution.vtu square_error_estimator_pcs_0_ts_4_t_1.000000.vtu displacement[1] displacement[1]
    #
    square_error_estimator_solution.vtu square_error_estimator_pcs_0_ts_4_t_1.000000.vtu sigma[0] sigma[0]
    square_error_estimator_solution.vtu square_error_estimator_pcs_0_ts_4_t_1.000000.vtu sigma[1] sigma[1]
    square_error_estimator_solution.vtu square_error_estimator_pcs_0_ts_4_t_1.000000.vtu sigma[2] sigma[2]
    square_error_estimator_solution.vtu square_error_estimator_pcs_0_ts_4_t_1.000000.vtu sigma[3] sigma[3]
    #
    square_error_estimator_solution.vtu square_error_estimator_pcs_0_ts_4_t_1.000000.vtu epsilon[0] epsilon[0]
    square_error_estimator_solution.vtu square_error_estimator_pcs_0_ts_4_t_1.000000.vtu epsilon[1] epsilon[1]
    square_error_estimator_solution.vtu square_error_estimator_pcs_0_ts_4_t_1.000000.vtu epsilon[2] epsilon[2]
    square_error_estimator_solution.vtu square_error_estimator_pcs_0_ts_4_t_1.000000.vtu epsilon[3] epsilon[3]
    )


# Simple test for the error estimator
# Stretching in y-direction with a non-uniform Dirichlet BC and local mesh refinement
# after global mesh refinement
# This test also verifies that the mesh refinement procedure comes to an end and does
# not run forever.
# Hint for verification: It is most important that the y-displacement at the top boundary
# matches the prescribed linear slope exactly!
AddTest(
    NAME DUNE_SmallDefLin_square_error_estimator_plus_global_refinement
    PATH Mechanics/LinearDUNE
    EXECUTABLE_ARGS square_error_estimator_plus_global_refinement.prj
    REQUIREMENTS NOT OGS_USE_MPI
    TESTER vtkdiff
    ABSTOL 1e-15 RELTOL 0
    DIFF_DATA
    square_error_estimator_plus_global_refinement_solution.vtu square_error_estimator_plus_global_refinement_pcs_0_ts_1_t_1.000000.vtu displacement[0] displacement[0]
    square_error_estimator_plus_global_refinement_solution.vtu square_error_estimator_plus_global_refinement_pcs_0_ts_1_t_1.000000.vtu displacement[1] displacement[1]
    #
    square_error_estimator_plus_global_refinement_solution.vtu square_error_estimator_plus_global_refinement_pcs_0_ts_1_t_1.000000.vtu sigma[0] sigma[0]
    square_error_estimator_plus_global_refinement_solution.vtu square_error_estimator_plus_global_refinement_pcs_0_ts_1_t_1.000000.vtu sigma[1] sigma[1]
    square_error_estimator_plus_global_refinement_solution.vtu square_error_estimator_plus_global_refinement_pcs_0_ts_1_t_1.000000.vtu sigma[2] sigma[2]
    square_error_estimator_plus_global_refinement_solution.vtu square_error_estimator_plus_global_refinement_pcs_0_ts_1_t_1.000000.vtu sigma[3] sigma[3]
    #
    square_error_estimator_plus_global_refinement_solution.vtu square_error_estimator_plus_global_refinement_pcs_0_ts_1_t_1.000000.vtu epsilon[0] epsilon[0]
    square_error_estimator_plus_global_refinement_solution.vtu square_error_estimator_plus_global_refinement_pcs_0_ts_1_t_1.000000.vtu epsilon[1] epsilon[1]
    square_error_estimator_plus_global_refinement_solution.vtu square_error_estimator_plus_global_refinement_pcs_0_ts_1_t_1.000000.vtu epsilon[2] epsilon[2]
    square_error_estimator_plus_global_refinement_solution.vtu square_error_estimator_plus_global_refinement_pcs_0_ts_1_t_1.000000.vtu epsilon[3] epsilon[3]
    )
