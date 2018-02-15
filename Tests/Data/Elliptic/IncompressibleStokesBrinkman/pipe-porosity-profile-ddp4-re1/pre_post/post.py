#!/usr/bin/env pvpython
# -*- coding: utf-8 -*-

from __future__ import print_function
import paraview.simple as ps
import cl_work.pre.via_paraview as pre

import sys
filename = sys.argv[1]

from pipe_params import *

center = (lx/2.0, ly/2.0, 0.0)
center_normal = (0.0, 1.0, 0.0)
print("center", center)
print("normal", center_normal)

# reader = ps.XMLUnstructuredGridReader(FileName=["out/pipe_pcs_0_ts_1_t_1.000000.vtu"])
# reader = ps.XMLUnstructuredGridReader(FileName=["out/pipe_pcs_0_ts_8_t_4.000000.vtu"])
# reader = ps.XMLUnstructuredGridReader(FileName=["out/pipe_pcs_0_ts_1_t_4.000000.vtu"])
reader = ps.XMLUnstructuredGridReader(FileName=[filename])

slice1 = pre.PlaneSlice(center, center_normal, Input=reader)

ps.SaveData("plot.csv", proxy=slice1, Precision=12) #, Scientific=True)
