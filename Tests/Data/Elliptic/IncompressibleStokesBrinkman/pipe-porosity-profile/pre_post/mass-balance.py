#!/usr/bin/env pvpython

from __future__ import print_function

import paraview
import numpy as np
import math

import sys

pvd_file = sys.argv[1]
velocity_field = sys.argv[2]

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# Create a new 'Render View'
view = GetActiveViewOrCreate('RenderView')

# create a new 'PVD Reader'
reader = PVDReader(FileName=pvd_file)

slice_ = Slice(Input=reader)
slice_.SliceType = 'Plane'
slice_.SliceType.Origin = (0.5, 3.5, 0.0)
slice_.SliceType.Normal = (0.0, 1.0, 0.0)

# create a new 'Calculator'
v_times_r = Calculator(Input=slice_)
v_times_r.ResultArrayName = 'v_times_r'
v_times_r.Function = '-2.0*3.14159265358979323846*coordsX * {}_Y'.format(velocity_field)

# create a new 'Integrate Variables'
integrateVariables1 = IntegrateVariables(Input=v_times_r)

Show(reader, view)
Show(slice_, view)

ts = reader.TimestepValues

print()
print("ideal:", math.pi/2.0)

for t in ts:
    print()
    view.ViewTime = t

    for y in np.arange(0.25, 5.0, 0.5):
        slice_.SliceType.Origin = (0.5, y, 0.0)
        Render()
        # input("x")

        res = paraview.servermanager.Fetch(integrateVariables1)

        data = res.GetPointData()
        total_flux = data.GetAbstractArray("v_times_r").GetComponent(0, 0)
        print(t, y, total_flux)
