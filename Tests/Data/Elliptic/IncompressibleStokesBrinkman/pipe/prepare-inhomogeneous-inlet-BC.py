#!/usr/bin/env pvpython
# -*- coding: utf-8 -*-

from __future__ import print_function

import inspect
import pprint

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--input", type=str, required=True)
parser.add_argument("--output", type=str, required=True)
parser.add_argument("--bc-field", type=str, required=True)
parser.add_argument("--total-flux", type=float, required=True)
parser.add_argument("--slice-normal", type=float, nargs=3, required=True)
parser.add_argument("--slice-origin", type=float, nargs=3, required=True)

args = parser.parse_args()

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'XML Unstructured Grid Reader'
xMLUnstructuredGridReader1 = XMLUnstructuredGridReader(FileName=[args.input])


def ProgrammableFilterWithScript(Script, parameters, *args, **kwargs):
    filt = ProgrammableFilter(*args, **kwargs)

    assert pprint.isreadable(parameters)
    ls = [ "parameters = " + pprint.pformat(parameters) + "\n\n" ]

    # properly indent script
    for i, l in enumerate(inspect.getsourcelines(Script)[0]):
        if i == 0:
            indent = len(l) - len(l.lstrip())
            ls.append("def script():\n")
        else:
            ls.append(l[indent:])

    ls.append("script()\n")

    filt.Script = "".join(ls)
    # print(filt.Script)
    filt.RequestInformationScript = ''
    filt.RequestUpdateExtentScript = ''
    filt.PythonPath = ''
    return filt


def EnumeratePoints(*args, **kwargs):
    def do_enumerate_points(self):
        data = self.GetInputDataObject(0, 0)

        # point_ids = vtk.vtkIdTypeArray()  ## TODO why does that not work?
        # point_ids = vtk.vtkLongArray()
        point_ids = vtk.vtkUnsignedLongArray()
        point_ids.SetName("OriginalSubsurfaceNodeIDs")
        point_ids.SetNumberOfComponents(1)
        N = data.GetPoints().GetNumberOfPoints()
        point_ids.SetNumberOfTuples(N)

        for i in range(N):
            point_ids.SetComponent(i, 0, i)

        self.GetOutputDataObject(0).GetPointData().AddArray(point_ids)

        ### Cell ids do not work like that:
        # cell_ids = vtk.vtkDoubleArray()
        # cell_ids.SetName("bulk_mesh_element_ids")
        # cell_ids.SetNumberOfComponents(1)
        # C = data.GetCellData().GetNumberOfTuples()
        # cell_ids.SetNumberOfTuples(N)

        # for i in range(C):
        #     cell_ids.SetComponent(i, 0, i)

        # self.GetOutputDataObject(0).GetCellData().AddArray(cell_ids)

    return ProgrammableFilterWithScript(do_enumerate_points, {}, *args, **kwargs)


def PlaneSlice(origin, normal, *args, **kwargs):
    slice_ = Slice(*args, **kwargs)
    slice_.SliceType = 'Plane'
    slice_.Triangulatetheslice = 0
    slice_.SliceOffsetValues = [0.0]

    # init the 'Plane' selected for 'SliceType'
    slice_.SliceType.Origin = origin
    slice_.SliceType.Normal = normal

    return slice_


def AsUnstructuredGrid(*args, **kwargs):
    # Trick from http://www.vtk.org/Wiki/VTK/Examples/Cxx/PolyData/PolyDataToUnstructuredGrid
    # for converting polydata to an unstructured grid
    return AppendDatasets(*args, **kwargs)


def ParabolicProfileInX_1Component(result, *args, **kwargs):
    def compute_profile():
        import numpy as np
        from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk

        data = self.GetInputDataObject(0, 0)

        xs = vtk_to_numpy(data.GetPoints().GetData())[:,0]
        profile = - (1.0 - xs**2)

        profile = numpy_to_vtk(profile, 1)
        profile.SetName(parameters["result"])

        self.GetOutputDataObject(0).GetPointData().AddArray(profile)

    parameters = {
            "result": result
            }
    return ProgrammableFilterWithScript(compute_profile, parameters, *args, **kwargs)


def ParabolicProfileInX_2Components(result, *args, **kwargs):
    def compute_profile():
        import numpy as np
        from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk

        data = self.GetInputDataObject(0, 0)

        xs = vtk_to_numpy(data.GetPoints().GetData())[:,0]
        profile = np.zeros((len(xs), 2))
        profile[:,1] = - (1.0 - xs**2)

        profile = numpy_to_vtk(profile, 1)
        profile.SetName(parameters["result"])

        self.GetOutputDataObject(0).GetPointData().AddArray(profile)

    parameters = {
            "result": result
            }
    return ProgrammableFilterWithScript(compute_profile, parameters, *args, **kwargs)



def LinearProfileInY_1Component(result, *args, **kwargs):
    def compute_profile():
        import numpy as np
        from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk

        data = self.GetInputDataObject(0, 0)

        ys = vtk_to_numpy(data.GetPoints().GetData())[:,1]
        profile = np.empty(len(ys))
        profile[:] = 4*ys

        profile = numpy_to_vtk(profile, 1)
        profile.SetName(parameters["result"])

        self.GetOutputDataObject(0).GetPointData().AddArray(profile)

    parameters = {
            "result": result
            }
    return ProgrammableFilterWithScript(compute_profile, parameters, *args, **kwargs)



enumerate_points = EnumeratePoints(Input=xMLUnstructuredGridReader1)
enumerate_points.CopyArrays = 0 # copy input to output arrays

slice1 = PlaneSlice(args.slice_origin, args.slice_normal, Input=enumerate_points)

appendDatasets1 = AsUnstructuredGrid(Input=slice1)

profile = ParabolicProfileInX_1Component(args.bc_field, Input=appendDatasets1)
profile.CopyArrays = 1 # copy input to output arrays

profile2 = ParabolicProfileInX_2Components(args.bc_field, Input=xMLUnstructuredGridReader1)
profile2.CopyArrays = 1 # copy input to output arrays

profile3 = LinearProfileInY_1Component("initial_pressure", Input=profile2)
profile3.CopyArrays = 1 # copy input to output arrays

# ----------------------------------------------------------------
# finally, restore active source
# SetActiveSource(profile)
# ----------------------------------------------------------------

# save data
SaveData(args.output, proxy=profile, DataMode='Binary',
    EncodeAppendedData=1,
    CompressorType='ZLib')

SaveData("ic_" + args.output, proxy=profile3, DataMode='Binary',
    EncodeAppendedData=1,
    CompressorType='ZLib')
