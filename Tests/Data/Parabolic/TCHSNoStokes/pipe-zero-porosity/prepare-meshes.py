#!/usr/bin/env pvpython
# -*- coding: utf-8 -*-

from __future__ import print_function
import paraview.simple as ps
import cl_work.pre.via_paraview as pre
import sys

OGSPATH = "/home/lehmannc/prog/ogs6/github-chleh-PRs/build-release-eigen/bin"

reader = ps.XMLUnstructuredGridReader(FileName=["pipe_quadratic.vtu"])

### Material Ids

def f_mat_ids(coords):
    import numpy as np
    return np.zeros(coords.shape[0]).astype("i")

def f_profile_T(coords):
    import numpy as np
    return 3e2 + 3e1 * np.sin(np.pi*coords[:,0]) * np.sin(np.pi*coords[:,1])

mat_ids = pre.CellFunction("MaterialIDs", f_mat_ids, Input=reader)
profile_T = pre.NodalFunction("temperature_0", f_profile_T, Input=mat_ids)
profile_T.CopyArrays = 1

# save data
ps.SaveData("pipe_ic.vtu", proxy=profile_T, DataMode='Binary',
    EncodeAppendedData=1,
    CompressorType='ZLib')

# reference solution

def f_profile_T_500(coords):
    import numpy as np
    return 3e2 + 3e1 * np.sin(np.pi*coords[:,0]) * np.sin(np.pi*coords[:,1]) \
            * np.exp(-2.0 * 0.00114443464762167888 * 500)

def f_profile_T_1000(coords):
    import numpy as np
    return 3e2 + 3e1 * np.sin(np.pi*coords[:,0]) * np.sin(np.pi*coords[:,1]) \
            * np.exp(-2.0 * 0.00114443464762167888 * 1000)

def f_p_ref(coords):
    import numpy as np
    return 1e5 * np.ones(coords.shape[0])

def f_xmV_ref(coords):
    import numpy as np
    return 0.5 * np.ones(coords.shape[0])

def f_v_Darcy_ref(coords):
    import numpy as np
    return np.zeros((coords.shape[0], 2))

T_ref = pre.NodalFunction("T_ref", f_profile_T_500, Input=reader)
p_ref = pre.NodalFunction("p_ref", f_p_ref, Input=T_ref)
p_ref.CopyArrays = 1
xmV_ref = pre.NodalFunction("xmV_ref", f_xmV_ref, Input=p_ref)
xmV_ref.CopyArrays = 1
v_ref = pre.NodalFunction("v_Darcy_ref", f_v_Darcy_ref, Input=xmV_ref)
v_ref.CopyArrays = 1

# save data
ps.SaveData("pipe_ref_500s.vtu", proxy=v_ref, DataMode='Binary',
    EncodeAppendedData=1,
    CompressorType='ZLib')


T_ref = pre.NodalFunction("T_ref", f_profile_T_1000, Input=reader)
p_ref = pre.NodalFunction("p_ref", f_p_ref, Input=T_ref)
p_ref.CopyArrays = 1
xmV_ref = pre.NodalFunction("xmV_ref", f_xmV_ref, Input=p_ref)
xmV_ref.CopyArrays = 1
v_ref = pre.NodalFunction("v_Darcy_ref", f_v_Darcy_ref, Input=xmV_ref)
v_ref.CopyArrays = 1

# save data
ps.SaveData("pipe_ref_1000s.vtu", proxy=v_ref, DataMode='Binary',
    EncodeAppendedData=1,
    CompressorType='ZLib')
