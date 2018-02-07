#!/usr/bin/env pvpython
# -*- coding: utf-8 -*-

from __future__ import print_function
import paraview.simple as ps
import cl_work.pre.via_paraview as pre

import os
import subprocess

OGSPATH = "/home/lehmannc/prog/ogs6/github-chleh-PRs/build-release-eigenlis/bin"

inlet_origin = (0.5, 1.0, 0.0)
inlet_normal = (0.0, 1.0, 0.0)

subprocess.check_call([os.path.join(OGSPATH, "generateStructuredMesh"),
    "-e", "quad",
    "--lx", "1",
    "--ly", "1",
    "--nx", "20",
    "--ny", "20",
    "-o", "pipe_linear.vtu"
    ])

subprocess.check_call([os.path.join(OGSPATH, "createQuadraticMesh"),
    "-i", "pipe_linear.vtu",
    "-o", "pipe.vtu"
    ])


reader = ps.XMLUnstructuredGridReader(FileName=["pipe.vtu"])

enumerate_points = pre.EnumeratePoints(Input=reader)
enumerate_points.CopyArrays = 0 # copy input to output arrays

slice1 = pre.PlaneSlice(inlet_origin, inlet_normal, Input=enumerate_points)

to_unstructured = pre.AsUnstructuredGrid(Input=slice1)

profile = pre.ParabolicProfileInX_1Component("velocity_y_inlet", Input=to_unstructured)
profile.CopyArrays = 1 # copy input to output arrays

# save data
ps.SaveData("pipe_bc_inlet.vtu", proxy=profile, DataMode='Binary',
    EncodeAppendedData=1,
    CompressorType='ZLib')
