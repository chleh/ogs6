AddTest(
    NAME TCHSStokes_pipe
    PATH Parabolic/TCHSStokes/pipe
    EXECUTABLE ogs
    EXECUTABLE_ARGS pipe.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    pipe_ref.vtu pipe_pcs_0_ts_1_t_1.000000.vtu v_ref darcy_velocity 1e-11 1e-15
    pipe_ref.vtu pipe_pcs_0_ts_1_t_1.000000.vtu p_ref pressure_interpolated 1e-15 1e-15
    pipe_ref.vtu pipe_pcs_0_ts_1_t_1.000000.vtu T_ref temperature_interpolated 1e-15 1e-6
    pipe_ref.vtu pipe_pcs_0_ts_1_t_1.000000.vtu x_ref mass_fraction_interpolated 1e-13 1e-13
)
