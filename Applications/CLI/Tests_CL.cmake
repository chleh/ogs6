### General Tests #######

# This test checks that
# * Dirichlet boundary conditions for temperature, pressure, mass fraction are set correctly
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
AddTest(
	NAME TESProcess-1D-0.1-10elem-by-comp-Dirichlet
	PATH Parabolic/TES_1D
	EXECUTABLE_ARGS line_0.1-10elem-by-comp-Dirichlet.prj
	WRAPPER ../../../../../ogs6-data-mine/ogs-wrapper.sh
	# TESTER vtkdiff
	# DIFF_DATA line_0.1_dirichlet_pcs_0_ts_1.vtu D1_left_bottom_N1_right Result
	DATA dummy
)


# This test checks that
# * Neumann BCs are put to the right component
AddTest(
	NAME TESProcess-2D-1x1-1e2elem-by-node-Neumann-heat-left
	PATH Parabolic/TES_2D
	EXECUTABLE_ARGS square_1x1-1e2elem-by-node-Neumann-heat-left.prj
	WRAPPER ../../../../../ogs6-data-mine/ogs-wrapper.sh
	# TESTER vtkdiff
	# DIFF_DATA line_0.1_dirichlet_pcs_0_ts_1.vtu D1_left_bottom_N1_right Result
	DATA dummy
)
AddTest(
	NAME TESProcess-2D-1x1-1e2elem-by-comp-Neumann-heat-left
	PATH Parabolic/TES_2D
	EXECUTABLE_ARGS square_1x1-1e2elem-by-comp-Neumann-heat-left.prj
	WRAPPER ../../../../../ogs6-data-mine/ogs-wrapper.sh
	# TESTER vtkdiff
	# DIFF_DATA line_0.1_dirichlet_pcs_0_ts_1.vtu D1_left_bottom_N1_right Result
	DATA dummy
)



### Physics Tests #######

AddTest(
	NAME TESProcess-1D-0.1-10elem-to-ads-equil
	PATH Parabolic/TES_1D
	EXECUTABLE_ARGS line_0.1-10elem-to-ads-equil.prj
	WRAPPER ../../../../../ogs6-data-mine/ogs-wrapper.sh
	# TESTER vtkdiff
	# DIFF_DATA line_0.1_dirichlet_pcs_0_ts_1.vtu D1_left_bottom_N1_right Result
	DATA dummy
)
