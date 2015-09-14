# This test checks that
# * dirichlet boundary conditions for temperature, pressure, mass fraction are set correctly
# * the process is correctly solved
# * the first and last timesteps are written to output files
AddTest(
	NAME TESProcess-1D-0.1-10elem-by-node-Dirichlet
	PATH Parabolic/TES_1D
	EXECUTABLE_ARGS line_0.1-10elem-by-node-Dirichlet.prj
	WRAPPER ../../../../../ogs6-data-mine/ogs-wrapper.sh
	# TESTER vtkdiff
	# DIFF_DATA line_0.1_dirichlet_pcs_0_ts_1.vtu D1_left_bottom_N1_right Result
	DATA dummy
)
