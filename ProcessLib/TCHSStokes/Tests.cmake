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
    pipe_ref.vtu pipe_pcs_0_ts_1_t_1.000000.vtu T_ref temperature_interpolated 1e-15 1e-14
    pipe_ref.vtu pipe_pcs_0_ts_1_t_1.000000.vtu x_ref mass_fraction_interpolated 1e-13 1e-13
)

AddTest(
    NAME TCHSStokes_pipe_porosity_profile_eff_viscosity_ddp10_re1
    PATH Parabolic/TCHSStokes/pipe-porosity-profile-eff-viscosity-ddp10-re1
    EXECUTABLE ogs
    EXECUTABLE_ARGS pipe.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    pipe_ref.vtu pipe_pcs_0_ts_1_t_1.000000.vtu pressure_interpolated pressure_interpolated 1e-15 1e-14
    pipe_ref.vtu pipe_pcs_0_ts_1_t_1.000000.vtu temperature_interpolated temperature_interpolated 1e-15 1e-13
    pipe_ref.vtu pipe_pcs_0_ts_1_t_1.000000.vtu mass_fraction_interpolated mass_fraction_interpolated 1e-15 1e-13
    pipe_ref.vtu pipe_pcs_0_ts_1_t_1.000000.vtu darcy_velocity darcy_velocity 1e-11 1e-15
)

AddTest(
    NAME TCHSStokes_pipe_zero_porosity
    PATH Parabolic/TCHSStokes/pipe-zero-porosity
    EXECUTABLE ogs
    EXECUTABLE_ARGS pipe.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    pipe_ref_500s.vtu pipe_pcs_0_ts_50_t_500.000000.vtu p_ref pressure_interpolated 1e-15 1e-15
    pipe_ref_500s.vtu pipe_pcs_0_ts_50_t_500.000000.vtu T_ref temperature_interpolated 1e-15 4e-4
    pipe_ref_500s.vtu pipe_pcs_0_ts_50_t_500.000000.vtu xmV_ref mass_fraction_interpolated 1e-15 1e-12
    pipe_ref_500s.vtu pipe_pcs_0_ts_50_t_500.000000.vtu v_Darcy_ref darcy_velocity 1e-15 1e-15
    #
    pipe_ref_1000s.vtu pipe_pcs_0_ts_100_t_1000.000000.vtu p_ref pressure_interpolated 1e-15 1e-15
    pipe_ref_1000s.vtu pipe_pcs_0_ts_100_t_1000.000000.vtu T_ref temperature_interpolated 1e-15 3e-4
    pipe_ref_1000s.vtu pipe_pcs_0_ts_100_t_1000.000000.vtu xmV_ref mass_fraction_interpolated 1e-15 1e-12
    pipe_ref_1000s.vtu pipe_pcs_0_ts_100_t_1000.000000.vtu v_Darcy_ref darcy_velocity 1e-15 1e-15
)

AddTest(
    NAME TCHSStokes_test_fluid_eos
    PATH Parabolic/TCHSStokes/test-fluid-eos
    EXECUTABLE ogs
    EXECUTABLE_ARGS pipe.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    pipe_ref.vtu pipe_pcs_0_ts_20_t_200.000000.vtu p_ref pressure_interpolated 1e-15 2e-7
    pipe_ref.vtu pipe_pcs_0_ts_20_t_200.000000.vtu T_ref temperature_interpolated 1e-15 1e-13
    pipe_ref.vtu pipe_pcs_0_ts_20_t_200.000000.vtu xmV_ref mass_fraction_interpolated 1e-15 1e-15
    pipe_ref.vtu pipe_pcs_0_ts_20_t_200.000000.vtu v_Darcy_ref darcy_velocity 1e-10 1e-15
)
