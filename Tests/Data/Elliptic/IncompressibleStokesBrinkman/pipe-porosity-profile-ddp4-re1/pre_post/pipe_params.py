pellet_diameter = 0.5

lx = 1.0
ly = 7.0
average_cell_size = 0.1

bed_start = 5.5
bed_end = 1.5
bed_radius = lx

refine_width_wall = min(lx, 1.5 * pellet_diameter) # 1.5 pellet diameters
refine_factor_wall = 50.0

refine_width_bed = 2*refine_width_wall
refine_factor_bed = 10.0

coord_tol = 1e-6

average_darcy_velocity = 0.5
