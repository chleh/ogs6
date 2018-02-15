#!/usr/bin/env pvpython

import paraview.simple as ps
import paraview
import sys
import inspect

fn = sys.argv[1]

reader = ps.XMLUnstructuredGridReader(FileName=[fn])

check_mesh = ps.ProgrammableFilter(Input=reader)
def do_check_mesh():
    import numpy as np
    from vtk.util.numpy_support import vtk_to_numpy
    data = self.GetInputDataObject(0, 0)
    coords = vtk_to_numpy(data.GetPoints().GetData())
    print(np.max(coords, axis=0))
check_mesh.Script = inspect.getsource(do_check_mesh) + "\n\ndo_check_mesh()"
check_mesh.CopyArrays = 1

paraview.servermanager.Fetch(check_mesh)
