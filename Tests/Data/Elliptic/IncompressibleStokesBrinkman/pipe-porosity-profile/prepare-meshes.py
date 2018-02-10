#!/usr/bin/env pvpython
# -*- coding: utf-8 -*-

from __future__ import print_function
import paraview.simple as ps
import cl_work.pre.via_paraview as pre

import os
import subprocess
from string import Template
import numpy as np

OGSPATH = "/home/lehmannc/prog/ogs6/github-chleh-PRs/build-release-eigenlis/bin"

lx = 1.0

ly = 5.0
bed_start = 4.0
bed_end = 1.0

nx = 80
ny = int(nx/lx*ly/2)
mx = 0.95

inlet_origin = (lx/2.0, ly, 0.0)
inlet_normal = (0.0, 1.0, 0.0)


gml = Template("""
<?xml version="1.0" encoding="ISO-8859-1"?>
<?xml-stylesheet type="text/xsl" href="OpenGeoSysGLI.xsl"?>

<OpenGeoSysGLI xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="http://www.opengeosys.org/images/xsd/OpenGeoSysGLI.xsd" xmlns:ogs="http://www.opengeosys.org">
    <name>pipe</name>
    <points>
        <point id="0" x="0" y="0" z="0"/>
        <point id="1" x="0" y="${ly}" z="0"/>
        <point id="2" x="${lx}" y="0" z="0"/>
        <point id="3" x="${lx}" y="${ly}" z="0"/>
    </points>

    <polylines>
        <polyline id="0" name="inner">
            <pnt>0</pnt>
            <pnt>1</pnt>
        </polyline>
        <polyline id="1" name="outer">
            <pnt>2</pnt>
            <pnt>3</pnt>
        </polyline>
        <polyline id="2" name="bottom">
            <pnt>0</pnt>
            <pnt>2</pnt>
        </polyline>
        <polyline id="3" name="top">
            <pnt>1</pnt>
            <pnt>3</pnt>
        </polyline>
    </polylines>
</OpenGeoSysGLI>
""").substitute(lx=lx, ly=ly)

with open("pipe.gml", "w") as fh:
    fh.write(gml)


### create mesh

subprocess.check_call([os.path.join(OGSPATH, "generateStructuredMesh"),
    "-e", "quad",
    "--lx", str(lx),
    "--ly", str(ly),
    "--nx", str(nx),
    "--ny", str(ny),
    "--mx", str(mx),
    "-o", "tmp_pipe_linear.vtu"
    ])

subprocess.check_call([os.path.join(OGSPATH, "createQuadraticMesh"),
    "-i", "tmp_pipe_linear.vtu",
    "-o", "tmp_pipe_quadratic_raw.vtu"
    ])


### create BCs

reader = ps.XMLUnstructuredGridReader(FileName=["tmp_pipe_quadratic_raw.vtu"])

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


### Material Ids

def f_mat_ids(coords):
    import numpy as np
    # True=1 if void, False=0 if bed
    void_cells = (coords[:,1] < bed_end) | (coords[:,1] > bed_start)
    void_val = int(0)
    bed_val = int(1)
    return np.choose(void_cells, (bed_val, void_val)).astype("i")

mat_ids = pre.CellFunction("MaterialIDs", f_mat_ids, Input=reader)

def f_profile_v(coords):
    import numpy as np
    profile = np.zeros((coords.shape[0], 2))
    profile[:,1] = - (1.0 - coords[:,0]**2)
    return profile

def f_profile_p(coords):
    return 4.0 * coords[:,1]

profile_v = pre.NodalFunction("initial_velocity", f_profile_v, Input=mat_ids)
profile_v.CopyArrays = 1
profile_p = pre.NodalFunction("initial_pressure", f_profile_p, Input=profile_v)
profile_p.CopyArrays = 1

# save data
ps.SaveData("pipe.vtu", proxy=profile_p, DataMode='Binary',
    EncodeAppendedData=1,
    CompressorType='ZLib')


xs = np.linspace(0, lx, 1001)
ys = f_profile_v_y(np.atleast_2d(xs).T)
A = 2 * np.pi * np.trapz(x=xs, y=xs)
Q = 2 * np.pi * np.trapz(x=xs, y=xs*ys)
print("cross sectional area:", A)
print("total flux:", Q, "m³/m²/s")
print("average velocity:", Q/A, "m/s")


### create reference solution

if False:
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


os.unlink("tmp_pipe_linear.vtu")
os.unlink("tmp_pipe_quadratic_raw.vtu")
