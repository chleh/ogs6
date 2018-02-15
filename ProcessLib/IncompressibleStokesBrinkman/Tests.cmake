# Incompressible Stokes Brinkman
AddTest(
    NAME IncompressibleStokesBrinkman_pipe
    PATH Elliptic/IncompressibleStokesBrinkman/pipe
    EXECUTABLE ogs
    EXECUTABLE_ARGS pipe.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    pipe_ref.vtu pipe_pcs_0_ts_1_t_1.000000.vtu v_ref darcy_velocity 1e-11 1e-15
    pipe_ref.vtu pipe_pcs_0_ts_1_t_1.000000.vtu p_ref pressure_interpolated 1e-9 1e-15
)


# Incompressible Stokes Brinkman Reynolds number 1
# Cf. Winterberg, M., Tsotsas, E., 2000. Impact of tube-to-particle-diameter ratio on pressure drop in packed beds. AIChE Journal 46, 1084–1088. doi:10.1002/aic.690460519
AddTest(
    NAME IncompressibleStokesBrinkman_pipe_poro_profile_Re1_ddp4
    PATH Elliptic/IncompressibleStokesBrinkman/pipe-porosity-profile-ddp4-re1
    EXECUTABLE ogs
    EXECUTABLE_ARGS pipe.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    pipe_ref.vtu pipe_pcs_0_ts_1_t_4.000000.vtu darcy_velocity darcy_velocity 1e-15 1e-15
    pipe_ref.vtu pipe_pcs_0_ts_1_t_4.000000.vtu pressure_interpolated pressure_interpolated 1e-15 1e-15
)


# Incompressible Stokes Brinkman Reynolds number 1000
# Cf. Winterberg, M., Tsotsas, E., 2000. Impact of tube-to-particle-diameter ratio on pressure drop in packed beds. AIChE Journal 46, 1084–1088. doi:10.1002/aic.690460519
AddTest(
    NAME IncompressibleStokesBrinkman_pipe_poro_profile_Re1000_ddp4
    PATH Elliptic/IncompressibleStokesBrinkman/pipe-porosity-profile-ddp4-re1000
    EXECUTABLE ogs
    EXECUTABLE_ARGS pipe.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    pipe_ref.vtu pipe_pcs_0_ts_1_t_4.000000.vtu darcy_velocity darcy_velocity 1e-15 1e-15
    pipe_ref.vtu pipe_pcs_0_ts_1_t_4.000000.vtu pressure_interpolated pressure_interpolated 1e-15 1e-15
)
