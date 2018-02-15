#!/usr/bin/env pvpython

import vtk
from vtk.util.numpy_support import vtk_to_numpy
import numpy as np

from pipe_params import *


def get_dndx(how, dndx_i, dndx_f, width):
    if how == "exp":
        b = 1/width * np.log(dndx_f/dndx_i)
        def dndx(x):
            return dndx_i * np.exp(b*x)
    elif how == "lin":
        b = (dndx_f - dndx_i) / width
        def dndx(x):
            return dndx_i + b*x
    return dndx



def get_coords(how, dndx_i, dndx_f, width):
    assert width > 0
    dndx = get_dndx(how, dndx_i, dndx_f, width)

    xe = 0
    xs = [ xe ]
    while xe < width:
        xe += dndx(xe)
        xs.append(xe)

    xs = np.array(xs)
    if xe > width:
        xs *= width/xe

    return xs

def patch(xs, x_i, x_f, how, dndx_i, dndx_f):
    width = x_f - x_i
    assert width > 0
    xs = np.hstack((xs[xs<x_i-coord_tol],
        get_coords(how, dndx_i, dndx_f, width) + x_i,
        xs[xs>x_f+coord_tol]))
    return xs

def patch_symm(xs, x_i, x_f, how, dndx_i, dndx_f):
    width = (x_f - x_i) / 2.0
    assert width > 0
    new_coords = get_coords(how, dndx_i, dndx_f, width)
    new_coords_rev = new_coords[-2::-1]
    xs = np.hstack((xs[xs<x_i-coord_tol],
        new_coords + x_i,
        2*width + x_i - new_coords_rev,
        xs[xs>x_f+coord_tol]))
    return xs


def create_mesh(file_name):
    xs = np.arange(0, lx + average_cell_size/2, average_cell_size)
    ys = np.arange(0, ly + average_cell_size/2, average_cell_size)


    # patch xs
    xs = patch(xs, lx - refine_width_wall, lx,
            "exp", average_cell_size, average_cell_size/refine_factor_wall)

    ys = patch_symm(ys, bed_start - refine_width_bed, bed_start + refine_width_bed,
            "lin", average_cell_size, average_cell_size/refine_factor_bed)

    ys = patch_symm(ys, bed_end - refine_width_bed, bed_end + refine_width_bed,
            "lin", average_cell_size, average_cell_size/refine_factor_bed)

    print("x range", min(xs), max(xs))
    print("y range", min(ys), max(ys))
    # print(np.diff(ys))

    assert all(np.diff(xs) > coord_tol)
    assert all(np.diff(ys) > coord_tol)

    create_quad_mesh(xs, ys, file_name)
    # create_tri_mesh(xs, ys, file_name)


def create_tri_mesh(xs, ys, file_name):
    nx = len(xs) - 1
    ny = len(ys) - 1

    points = vtk.vtkPoints()
    points.SetDataTypeToDouble()
    i = 0
    for y in ys:
        for x in xs:
            points.InsertPoint(i, (x, y, 0))
            i += 1

    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.Allocate(2 * nx * ny)

    for iy in xrange(ny):
        for ix in xrange(nx):
            v0 = iy*(nx + 1) + ix
            v1 = v0 + 1
            v2 = v0 + nx+1
            v3 = v1 + nx+1

            ugrid.InsertNextCell(vtk.VTK_TRIANGLE, 3, (v0, v1, v3))
            ugrid.InsertNextCell(vtk.VTK_TRIANGLE, 3, (v0, v3, v2))

    ugrid.SetPoints(points)

    writer = vtk.vtkXMLUnstructuredGridWriter()
    # writer.SetDataMode(vtk.vtkXMLWriter.Ascii)
    writer.SetFileName(file_name)
    writer.SetInputDataObject(ugrid)
    writer.Write()

def create_quad_mesh(xs, ys, file_name):
    # Create a rectilinear grid by defining three arrays specifying the
    # coordinates in the x-y-z directions.
    xCoords = vtk.vtkDoubleArray()
    for i in xs:
        xCoords.InsertNextValue(i)

    yCoords = vtk.vtkDoubleArray()
    for i in ys:
        yCoords.InsertNextValue(i)

    zCoords = vtk.vtkDoubleArray()
    zCoords.InsertNextValue(0)

    # The coordinates are assigned to the rectilinear grid. Make sure that
    # the number of values in each of the XCoordinates, YCoordinates,
    # and ZCoordinates is equal to what is defined in SetDimensions().
    #
    rgrid = vtk.vtkRectilinearGrid()
    rgrid.SetDimensions(len(xs), len(ys), 1)
    rgrid.SetXCoordinates(xCoords)
    rgrid.SetYCoordinates(yCoords)
    rgrid.SetZCoordinates(zCoords)

    print("last point in rgrid", rgrid.GetPoint(len(xs)*len(ys)-1))

    if True:
        # convert to unstructured grid
        f_append = vtk.vtkAppendFilter()
        f_append.SetOutputPointsPrecision(1) # double precision
        f_append.SetInputDataObject(rgrid)

        f_append.Update()
        print(f_append.GetOutput().GetPoints())
        print("last point in appended grid:",
                f_append.GetOutput().GetPoints().GetPoint(len(xs) * len(ys) - 1))
        points = vtk_to_numpy(f_append.GetOutput().GetPoints().GetData())
        print(np.min(points, axis=0), np.max(points, axis=0))

        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetFileName(file_name)
        writer.SetInputConnection(f_append.GetOutputPort())
        writer.Write()
    else:
        writer = vtk.vtkXMLRectilinearGridWriter()
        writer.SetFileName(file_name)
        writer.SetInputDataObject(rgrid)
        writer.Write()
