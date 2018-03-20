AddTest(
    NAME TCHSNoStokes_pipe_zero_porosity
    PATH Parabolic/TCHSNoStokes/pipe-zero-porosity
    EXECUTABLE ogs
    EXECUTABLE_ARGS pipe.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    pipe_ref_500s.vtu pipe_pcs_0_ts_50_t_500.000000.vtu p_ref pressure 1e-15 1e-12
    pipe_ref_500s.vtu pipe_pcs_0_ts_50_t_500.000000.vtu T_ref temperature 1e-15 5e-4
    pipe_ref_500s.vtu pipe_pcs_0_ts_50_t_500.000000.vtu xmV_ref mass_fraction 1e-15 1e-12
    #
    pipe_ref_1000s.vtu pipe_pcs_0_ts_100_t_1000.000000.vtu p_ref pressure 1e-15 1e-12
    pipe_ref_1000s.vtu pipe_pcs_0_ts_100_t_1000.000000.vtu T_ref temperature 1e-15 3e-4
    pipe_ref_1000s.vtu pipe_pcs_0_ts_100_t_1000.000000.vtu xmV_ref mass_fraction 1e-15 1e-12
)

AddTest(
    NAME TCHSNoStokes_test_fluid_eos
    PATH Parabolic/TCHSNoStokes/test-fluid-eos
    EXECUTABLE ogs
    EXECUTABLE_ARGS pipe.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    pipe_ref.vtu pipe_pcs_0_ts_20_t_200.000000.vtu p_ref pressure 1e-15 2e-7
    pipe_ref.vtu pipe_pcs_0_ts_20_t_200.000000.vtu T_ref temperature 1e-15 1e-12
    pipe_ref.vtu pipe_pcs_0_ts_20_t_200.000000.vtu xmV_ref mass_fraction 1e-15 1e-14
)

AddTest(
    NAME TCHSNoStokes_reaction_sealed_container
    PATH Parabolic/TCHSNoStokes/reaction-sealed-container
    EXECUTABLE ogs
    EXECUTABLE_ARGS pipe.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    pipe_ref_0.2s.vtu pipe_pcs_0_ts_200_t_0.200000.vtu p_ref pressure        1e-15 6e-5
    pipe_ref_0.2s.vtu pipe_pcs_0_ts_200_t_0.200000.vtu T_ref temperature     1e-15 3e-5
    pipe_ref_0.2s.vtu pipe_pcs_0_ts_200_t_0.200000.vtu xmV_ref mass_fraction 1e-15 4e-5
    pipe_ref_0.2s.vtu pipe_pcs_0_ts_200_t_0.200000.vtu solid_density_ref solid_density    1e-15 1e-14
    pipe_ref_0.2s.vtu pipe_pcs_0_ts_200_t_0.200000.vtu reaction_rate_ref reaction_rate    1e-15 1e-14
    #
    pipe_ref_0.4s.vtu pipe_pcs_0_ts_400_t_0.400000.vtu p_ref pressure        1e-15 2e-4
    pipe_ref_0.4s.vtu pipe_pcs_0_ts_400_t_0.400000.vtu T_ref temperature     1e-15 5e-5
    pipe_ref_0.4s.vtu pipe_pcs_0_ts_400_t_0.400000.vtu xmV_ref mass_fraction 1e-15 9e-5
    pipe_ref_0.4s.vtu pipe_pcs_0_ts_400_t_0.400000.vtu solid_density_ref solid_density    1e-15 1e-14
    pipe_ref_0.4s.vtu pipe_pcs_0_ts_400_t_0.400000.vtu reaction_rate_ref reaction_rate    1e-15 1e-14
    #
    pipe_ref_0.6s.vtu pipe_pcs_0_ts_600_t_0.600000.vtu p_ref pressure        1e-15 2e-4
    pipe_ref_0.6s.vtu pipe_pcs_0_ts_600_t_0.600000.vtu T_ref temperature     1e-15 7e-5
    pipe_ref_0.6s.vtu pipe_pcs_0_ts_600_t_0.600000.vtu xmV_ref mass_fraction 1e-15 2e-4
    pipe_ref_0.6s.vtu pipe_pcs_0_ts_600_t_0.600000.vtu solid_density_ref solid_density    1e-15 1e-14
    pipe_ref_0.6s.vtu pipe_pcs_0_ts_600_t_0.600000.vtu reaction_rate_ref reaction_rate    1e-15 1e-14
    #
    pipe_ref_0.8s.vtu pipe_pcs_0_ts_800_t_0.800000.vtu p_ref pressure        1e-15 3e-4
    pipe_ref_0.8s.vtu pipe_pcs_0_ts_800_t_0.800000.vtu T_ref temperature     1e-15 8e-5
    pipe_ref_0.8s.vtu pipe_pcs_0_ts_800_t_0.800000.vtu xmV_ref mass_fraction 1e-15 2e-4
    pipe_ref_0.8s.vtu pipe_pcs_0_ts_800_t_0.800000.vtu solid_density_ref solid_density    1e-15 1e-13
    pipe_ref_0.8s.vtu pipe_pcs_0_ts_800_t_0.800000.vtu reaction_rate_ref reaction_rate    1e-15 1e-13
    #
    pipe_ref_1.0s.vtu pipe_pcs_0_ts_1000_t_1.000000.vtu p_ref pressure        1e-15 4e-4
    pipe_ref_1.0s.vtu pipe_pcs_0_ts_1000_t_1.000000.vtu T_ref temperature     1e-15 1e-4
    pipe_ref_1.0s.vtu pipe_pcs_0_ts_1000_t_1.000000.vtu xmV_ref mass_fraction 1e-15 3e-4
    pipe_ref_1.0s.vtu pipe_pcs_0_ts_1000_t_1.000000.vtu solid_density_ref solid_density    1e-15 1e-13
    pipe_ref_1.0s.vtu pipe_pcs_0_ts_1000_t_1.000000.vtu reaction_rate_ref reaction_rate    1e-15 1e-13
    #
)
