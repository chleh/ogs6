#!/usr/bin/env pvpython

from __future__ import print_function

import sys
from paraview.simple import *
import inspect
import pprint

fn1 = sys.argv[1]
field1 = sys.argv[2]
fn2 = sys.argv[3]
field2 = sys.argv[4]
out_prefix = sys.argv[5]

reader1 = XMLUnstructuredGridReader(FileName = [ fn1 ])
reader2 = XMLUnstructuredGridReader(FileName = [ fn2 ])


resampleWithDataset = ResampleWithDataset(Input=reader2,
    Source=reader1)


diff_filter = ProgrammableFilter(Input = [ reader1, resampleWithDataset ])
diff_filter.SetPropertyWithName("Parameters", [
    "field1", pprint.pformat(field1),
    "field2", pprint.pformat(field2),
    "out_prefix", pprint.pformat(out_prefix)
    ])

def do_diff():
    import numpy as np
    from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
    #import matplotlib.pyplot as plt

    assert self.GetNumberOfInputPorts() == 1
    assert self.GetNumberOfInputConnections(0) == 2
    print(field1)
    print(field2)

    data1 = self.GetInputDataObject(0, 0)
    data2 = self.GetInputDataObject(0, 1)

    pd1 = data1.GetPointData()
    pd2 = data2.GetPointData()

    f1_vtk = pd1.GetAbstractArray(field1)
    f2_vtk = pd2.GetAbstractArray(field2)

    f1 = vtk_to_numpy(f1_vtk)
    f2 = vtk_to_numpy(f2_vtk)

    if f1_vtk.GetNumberOfComponents() == 1:
        f1 = np.atleast_2d(f1).T
    if f2_vtk.GetNumberOfComponents() == 1:
        f2 = np.atleast_2d(f2).T

    f1_pts, f1_comp = f1.shape
    f2_pts, f2_comp = f2.shape
    assert f1_pts == f2_pts

    if f1_comp < f2_comp:
        f1 = np.hstack((f1, np.zeros((f1_pts, f2_comp-f1_comp))))
    if f2_comp < f1_comp:
        f2 = np.hstack((f2, np.zeros((f2_pts, f1_comp-f2_comp))))

    abs_diff = f1 - f2
    magnitude = 0.5 * (np.abs(f1) + np.abs(f2)) + 1e-16

    if False:
        if f1_comp == f2_comp:
            abs_diff = f1 - f2
            magnitude = 0.5 * (np.abs(f1) + np.abs(f2))
        elif f1_comp > f2_comp:
            abs_diff = f1
            abs_diff[:, 0:f2_comp] -= f2

            magnitude = np.empty_like(f1)
            magnitude[:, 0:f2_comp] = 0.5 * (np.abs(f1[:, 0:f2_comp]) + np.abs(f2))
            magnitude[:, f2_comp:f1_comp] = np.abs(f1[:, f2_comp:f1_comp])
        else:
            abs_diff = -f2
            abs_diff[:, 0:f1_comp] += f1

            magnitude = np.empty_like(f2)
            magnitude[:, 0:f1_comp] = 0.5 * (np.abs(f2[:, 0:f1_comp]) + np.abs(f1))
            magnitude[:, f1_comp:f2_comp] = np.abs(f2[:, f1_comp:f2_comp])

    rel_diff = abs_diff / magnitude

    points = vtk_to_numpy(data1.GetPoints().GetData())

    # print("abs diff min", np.min(abs_diff, axis=0), "max", np.max(abs_diff, axis=0))
    # print("abs diff argmin", np.argmin(abs_diff, axis=0), "argmax", np.argmax(abs_diff, axis=0))

    print("\nvalues, indices and positions of minimum absolute error:")
    for d, i in enumerate(np.argmin(abs_diff, axis=0)):
        print("  d={} min={}\t{}\t{}".format(d, abs_diff[i,d], i, points[i,:]))
        # print("  {}\t{}\t{}".format(abs_diff[i,d], points[i,:]))

    print("\nvalues, indices and positions of maximum absolute error:")
    for d, i in enumerate(np.argmax(abs_diff, axis=0)):
        print("  d={} max={}\t{}\t{}".format(d, abs_diff[i,d], i, points[i,:]))

    # print("rel diff min", np.min(rel_diff, axis=0), "max", np.max(rel_diff, axis=0))
    # print("rel diff argmin", np.argmin(rel_diff, axis=0), "argmax", np.argmax(rel_diff, axis=0))

    print("\nvalues, indices and positions of minimum relative error:")
    for d, i in enumerate(np.argmin(rel_diff, axis=0)):
        print("  d={} min={}\t{}\t{}".format(d, rel_diff[i,d], i, points[i,:]))

    print("\nvalues, indices and positions of maximum relative error:")
    for d, i in enumerate(np.argmax(rel_diff, axis=0)):
        print("  d={} max={}\t{}\t{}".format(d, rel_diff[i,d], i, points[i,:]))

    data_range = np.max(abs_diff, axis=0) - np.min(abs_diff, axis=0)
    for d in xrange(data_range.size):
        if data_range[d] == 0: continue
        hist, edges = np.histogram(abs_diff[:, d], bins=int(np.sqrt(f1_pts))+1)
        np.savez_compressed(out_prefix + "_diff_hist_abs_diff_comp_{}.npz".format(d), hist=hist, edges=edges)
        # fig, ax = plt.subplots()
        # ax.hist(abs_diff[:,d], bins=int(np.sqrt(f1_pts))+1)
        # ax.set_ylim(0.5)
        # ax.set_yscale("log")
        # fig.savefig("diff_abs_diff_comp_{}.png".format(d))
        # plt.close(fig)

    data_range = np.max(rel_diff, axis=0) - np.min(rel_diff, axis=0)
    for d in xrange(data_range.size):
        if data_range[d] == 0: continue
        hist, edges = np.histogram(rel_diff[:, d], bins=int(np.sqrt(f1_pts))+1)
        np.savez_compressed(out_prefix + "_diff_hist_rel_diff_comp_{}.npz".format(d), hist=hist, edges=edges)


    abs_diff = numpy_to_vtk(abs_diff, 1)
    rel_diff = numpy_to_vtk(rel_diff, 1)

    abs_diff.SetName("abs_diff_F1_{}_minus_F2_{}".format(field1, field2))
    rel_diff.SetName("rel_diff_F1_{}_minus_F2_{}".format(field1, field2))

    f1_vtk.SetName("F1_"+field1)
    f2_vtk.SetName("F2_"+field2)

    self.GetOutputDataObject(0).GetPointData().AddArray(abs_diff)
    self.GetOutputDataObject(0).GetPointData().AddArray(rel_diff)
    self.GetOutputDataObject(0).GetPointData().AddArray(f1_vtk)
    self.GetOutputDataObject(0).GetPointData().AddArray(f2_vtk)
    self.GetOutputDataObject(0).GetPointData().AddArray(pd2.GetAbstractArray("vtkValidPointMask"))





diff_filter.Script = inspect.getsource(do_diff) + "\n\ndo_diff()"
diff_filter.CopyArrays = 0

SaveData(out_prefix + "_diff.vtu", proxy=diff_filter, DataMode='Binary',
    EncodeAppendedData=1,
    CompressorType='ZLib')
