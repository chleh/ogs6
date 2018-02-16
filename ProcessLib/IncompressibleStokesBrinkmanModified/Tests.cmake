# Incompressible Stokes Brinkman
AddTest(
    NAME IncompressibleStokesBrinkmanModified_pipe
    PATH Elliptic/IncompressibleStokesBrinkmanModified/pipe
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
# VDI Waermeatlas
AddTest(
    NAME IncompressibleStokesBrinkmanModified_pipe_poro_profile_eff_viscosity_Re1_ddp10
    PATH Elliptic/IncompressibleStokesBrinkmanModified/pipe-porosity-profile-eff-viscosity-ddp10-re1
    EXECUTABLE ogs
    EXECUTABLE_ARGS pipe.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    pipe_ref.vtu pipe_pcs_0_ts_1_t_4.000000.vtu darcy_velocity darcy_velocity 1e-12 1e-15
    pipe_ref.vtu pipe_pcs_0_ts_1_t_4.000000.vtu pressure_interpolated pressure_interpolated 1e-10 1e-15
)
