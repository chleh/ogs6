#!/usr/bin/env pvpython
# -*- coding: utf-8 -*-

from __future__ import print_function
import paraview.simple as ps
import cl_work.pre.via_paraview as pre

import os
import subprocess
from string import Template
import numpy as np

bed_start = 5.5
bed_end = 1.5

### create mesh

reader = ps.XMLUnstructuredGridReader(FileName=["pipe_orig.vtu"])

def f_poro(coords):
    import numpy as np
    # True=1 if void, False=0 if bed
    void_cells = (coords[:,1] < bed_end) | (coords[:,1] > bed_start)
    void_val = 1.0
    bed_val = 0.5
    return np.choose(void_cells, (bed_val, void_val))

poro = pre.CellFunction("porosity", f_poro, Input=reader)

def f_mat(coords):
    import numpy as np
    return np.zeros(coords.shape[0]).astype("i")

mat = pre.CellFunction("MaterialIDs", f_mat, Input=poro)
mat.CopyArrays = 1

# save data
ps.SaveData("pipe_poro.vtu", proxy=mat, DataMode='Binary',
    EncodeAppendedData=1,
    CompressorType='ZLib')
