AddTest(
    NAME TCHSStokes_test1
    PATH Parabolic/TCHSStokes/test1
    EXECUTABLE ogs
    EXECUTABLE_ARGS pipe.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    pipe_ref.vtu pipe_pcs_0_ts_1_t_4.000000.vtu darcy_velocity darcy_velocity 1e-12 1e-15
    pipe_ref.vtu pipe_pcs_0_ts_1_t_4.000000.vtu pressure_interpolated pressure_interpolated 1e-10 1e-15
)
