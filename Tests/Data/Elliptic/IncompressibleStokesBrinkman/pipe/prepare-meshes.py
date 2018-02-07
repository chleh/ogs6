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


### create mesh

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


### create BCs

reader = ps.XMLUnstructuredGridReader(FileName=["pipe.vtu"])

enumerate_points = pre.EnumeratePoints(Input=reader)
enumerate_points.CopyArrays = 0 # copy input to output arrays

slice1 = pre.PlaneSlice(inlet_origin, inlet_normal, Input=enumerate_points)

to_unstructured = pre.AsUnstructuredGrid(Input=slice1)

def f_profile_v_y(coords):
    return - (1.0 - coords[:,0]**2)

profile = pre.NodalFunction("velocity_y_inlet", f_profile_v_y, Input=to_unstructured)
profile.CopyArrays = 1 # copy input to output arrays

# save data
ps.SaveData("pipe_bc_inlet.vtu", proxy=profile, DataMode='Binary',
    EncodeAppendedData=1,
    CompressorType='ZLib')


### create reference solution

def f_profile_v(coords):
    import numpy as np
    profile = np.zeros((coords.shape[0], 2))
    profile[:,1] = - (1.0 - coords[:,0]**2)
    return profile

def f_profile_p(coords):
    return 4.0 * coords[:,1]

profile_v = pre.NodalFunction("v_ref", f_profile_v, Input=reader)
profile_p = pre.NodalFunction("p_ref", f_profile_p, Input=profile_v)
profile_p.CopyArrays = 1

ps.SaveData("pipe_ref.vtu", proxy=profile_p, DataMode='Binary',
    EncodeAppendedData=1,
    CompressorType='ZLib')
