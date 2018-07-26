#!/usr/bin/vtkpython

from __future__ import print_function
import vtk
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import numpy as np

import matplotlib.pyplot as plt

pvd_file = "out/hertz_pcs_0.pvd"


R12 = 1.0
R = R12 / 2.0
nu_12 = 0.3
E_12 = 1.0

lambda_ = E_12 * nu_12 / (1+nu_12) / (1-2*nu_12)
mu = E_12 / 2.0 / (1+nu_12)
G_12 = E_12 / 2.0 / (1 + nu_12)

kappa = 0.5 * G_12 / R / (1-nu_12)
print("kappa:", kappa)

C = lambda_ * np.matrix([
    1, 1, 1, 0,
    1, 1, 1, 0,
    1, 1, 1, 0,
    0, 0, 0, 0 ]).reshape(4, 4) \
            + 2 * mu * np.identity(4)

def p_contact(r, r_contact):
    return kappa * np.sqrt(r_contact**2 - r**2)


### helpers ##############################################

import os
try:
    import xml.etree.cElementTree as ET
except:
    import xml.etree.ElementTree as ET


def relpathfrom(origin, relpath):
    if os.path.isabs(relpath):
        return relpath
    return os.path.join(origin, relpath)

def read_pvd_file(fn):
    try:
        path = fn.name
    except AttributeError:
        path = fn
    pathroot = os.path.dirname(path)
    pvdtree = ET.parse(fn)
    node = pvdtree.getroot()
    if node.tag != "VTKFile": return None, None
    children = list(node)
    if len(children) != 1: return None, None
    node = children[0]
    if node.tag != "Collection": return None, None

    ts = []
    fs = []

    for child in node:
        if child.tag != "DataSet": return None, None
        ts.append(float(child.get("timestep")))
        fs.append(relpathfrom(pathroot, child.get("file")))

    return ts, fs

### helpers end ##########################################


def get_y_top(t):
    return 1.0 - 0.005 * t

ts, fns = read_pvd_file(pvd_file)

reader = vtk.vtkXMLUnstructuredGridReader()

strainFilter = vtk.vtkCellDerivatives()
strainFilter.SetVectorModeToPassVectors()
strainFilter.SetTensorModeToComputeStrain()

warpVector = vtk.vtkWarpVector()
warpVector.SetInputConnection(strainFilter.GetOutputPort())

destroyTopology = vtk.vtkShrinkFilter()
destroyTopology.SetShrinkFactor(1.0)
# destroyTopology.SetInputConnection(strainFilter.GetOutputPort())
destroyTopology.SetInputConnection(reader.GetOutputPort())

# cell2point = vtk.vtkCellDataToPointData()
# cell2point.SetInputConnection(destroyTopology.GetOutputPort())
# cell2point.PassCellDataOff()

plane = vtk.vtkPlane()
# plane.SetOrigin(0, 0, 0)
plane.SetNormal(0, 1, 0)

cutter = vtk.vtkCutter()
cutter.SetCutFunction(plane)
cutter.SetInputConnection(warpVector.GetOutputPort())

writer = vtk.vtkXMLUnstructuredGridWriter()

ws = []
rs_contact = []

fig, ax = plt.subplots()

for t, fn in zip(ts, fns):
    # print(t, fn)
    reader.SetFileName(fn)
    reader.Update()
    grid = reader.GetOutput()

    disp_2d = vtk_to_numpy(grid.GetPointData().GetArray("displacement"))
    disp_3d = np.zeros((disp_2d.shape[0], 3))
    disp_3d[:,(0,1)] = disp_2d
    disp_3d_vtk = numpy_to_vtk(disp_3d, 1)
    disp_3d_vtk.SetName("u")

    grid.GetPointData().AddArray(disp_3d_vtk)
    grid.GetPointData().SetActiveVectors("u")

    if True:
        # compute strain
        def strain_triangle_axi(cell, point_data, strain_data):
            cell_pts = np.matrix(vtk_to_numpy(cell.GetPoints().GetData())[:,:-1])
            assert cell_pts.shape[0] == 3  # triangles
            assert point_data.shape[1] == 2  # 2D vector field
            node_ids = [ cell.GetPointId(i) for i in range(cell.GetNumberOfPoints()) ]

            T = np.matrix(np.empty((2,2)))
            T[:,0] = (cell_pts[0,:] - cell_pts[2,:]).T
            T[:,1] = (cell_pts[1,:] - cell_pts[2,:]).T
            T_inv = np.linalg.inv(T)

            dl1 = T_inv[0,:]  # row 0
            dl2 = T_inv[1,:]  # row 1

            for node in range(3):
                l1, l2 = T_inv * (cell_pts[node, :].T - cell_pts[2, :].T)
                assert -1e-15 < l1 and 1 + 1e-15 > l1 \
                        and -1e-15 < l2 and 1+ 1e-15 > l2

            grad = np.empty((2,2))
            for comp in range(2):
                nodal_values = point_data[node_ids, comp]
                # nodal_values = cell_pts[:, comp].flat
                # if t > 0 and cell_pts[0,1] > 0.95 and comp == 1:
                #     print(nodal_values[0])
                grad[comp,:] = dl1 * nodal_values[0] \
                        + dl2 * nodal_values[1] \
                        - (dl1 + dl2) * nodal_values[2]

            # if t > 0 and cell_pts[0,1] > 0.95:
            #     print(grad)

            strain = 0.5 * (grad + grad.T)  # rr, rz, zr,zz components

            for node in range(3):
                r = cell_pts[node, 0]
                node_id = node_ids[node]

                if r == 0:
                    dvdr = grad[0,0]
                    v_over_r = dvdr
                else:
                    v_over_r = point_data[node_id, 0] / r

                strain_kelvin = np.array([
                    strain[0,0], strain[1,1], v_over_r,
                    strain[0,1] * np.sqrt(2.0)
                    ])
                strain_data[node_id, :] = strain_kelvin

        def computeStrain():
            destroyTopology.Update()
            grid = destroyTopology.GetOutput()

            disp_2d = vtk_to_numpy(grid.GetPointData().GetArray("displacement"))
            strain_kelvin = np.empty((disp_2d.shape[0], 4))

            n_cells = grid.GetNumberOfCells()
            for c in xrange(n_cells):
                cell = grid.GetCell(c)
                strain_triangle_axi(cell, disp_2d, strain_kelvin)

            strain_kelvin_vtk = numpy_to_vtk(strain_kelvin, 1)
            strain_kelvin_vtk.SetName("strain_post_kelvin")
            grid.GetPointData().AddArray(strain_kelvin_vtk)

            strain = strain_kelvin.copy()
            strain[:,3] /= np.sqrt(2.0)
            strain_vtk = numpy_to_vtk(strain, 1)
            strain_vtk.SetName("strain_post")
            grid.GetPointData().AddArray(strain_vtk)

            writer.SetInputData(grid)
            writer.SetFileName(os.path.join(
                os.path.dirname(fn), "post_{:.0f}.vtu".format(t)))
            writer.Write()

        computeStrain()

    strainFilter.SetInputData(grid)

    # compute stress for debugging
    if False:
        strainFilter.Update()
        grid = strainFilter.GetOutput()

        # FIXME problematic due to rotational symmetry
        strain = vtk_to_numpy(grid.GetCellData().GetArray("Strain"))
        strain_kelv = strain[:,(0, 4, 8, 1)]
        strain_kelv[:, 3] *= np.sqrt(2.0)

        stress_kelv = np.empty_like(strain_kelv)
        for c, eps in enumerate(strain_kelv):
            # print(C * np.atleast_2d(eps).T)
            stress_kelv[c, :] = (C * np.atleast_2d(eps).T).flat

        stress_symm_tensor = stress_kelv.copy()
        stress_symm_tensor[:,3] /= np.sqrt(2.0)

        stress_symm_tensor_vtk = numpy_to_vtk(stress_symm_tensor, 1)
        stress_symm_tensor_vtk.SetName("stress_post")

        grid.GetCellData().AddArray(stress_symm_tensor_vtk)

        writer.SetInputData(grid)
        writer.SetFileName(os.path.join(
            os.path.dirname(fn), "post_{:.0f}.vtu".format(t)))
        writer.Write()

    # warpVector.SetInputData(grid)
    warpVector.Update()
    grid = warpVector.GetOutput()

    xmin, xmax, ymin, ymax, zmin, zmax = grid.GetBounds()
    y_top = get_y_top(t)
    assert ymax == y_top
    ws.append(2 * (1.0 - y_top))

    # determine top boundary
    coords = vtk_to_numpy(grid.GetPoints().GetData())
    idcs_top_boundary = np.where(coords[:,1] > y_top - 1e-8)[0]
    xs_top_boundary = coords[idcs_top_boundary, 0]
    r_contact = max(xs_top_boundary) - min(xs_top_boundary)
    print("radius of contact area:", r_contact)

    rs_contact.append(r_contact)

    plane.SetOrigin(0, y_top, 0)
    cutter.Update()
    grid = cutter.GetOutput()

    # convert 3x3 strain tensor to 4x1 Kelvin vector
    strain = vtk_to_numpy(grid.GetCellData().GetArray("Strain"))
    strain_kelv = strain[:,(0, 4, 8, 1)]
    strain_kelv[:, 3] *= np.sqrt(2)

    stress_kelv = np.empty_like(strain_kelv)
    for c, eps in enumerate(strain_kelv):
        # print(C * np.atleast_2d(eps).T)
        stress_kelv[c, :] = (C * np.atleast_2d(eps).T).flat

    # print(stress_kelv)

    # TODO cell centers
    xs = vtk_to_numpy(grid.GetPoints().GetData())[:-1, 0]
    idcs = np.argsort(xs)
    # h, = ax.plot(xs[idcs], stress_kelv[:, 1][idcs])  # sigma_yy
    h, = ax.step(xs[idcs], stress_kelv[:, 1][idcs])  # sigma_yy
    assert abs(max(xs) - r_contact) < 1e-8
    # assert abs(min(xs)) < 1e-8

    rs = np.linspace(0, r_contact, 200)
    ax.plot(rs, -p_contact(rs, r_contact), color=h.get_color(), ls=":")

    rs = vtk_to_numpy(grid.GetPoints().GetData())[:,0]
    stress_yy = vtk_to_numpy(grid.GetPointData().GetArray("sigma"))[:,1]
    idcs = np.argsort(rs)
    ax.plot(rs[idcs], stress_yy[idcs], color=h.get_color(), ls="--")

fig.savefig("stress_at_contact.png")
plt.close(fig)


fig, ax = plt.subplots()
ax.scatter(ws, rs_contact)

ws_ref = np.linspace(0, max(ws), 200)
rs_ref = np.sqrt(ws_ref * R)

ax.plot(ws_ref, rs_ref)

ax.set_xlabel("displacement of sphere centers / m")
ax.set_ylabel("radius of contact areas / m")

fig.savefig("contact_radii.png")
plt.close(fig)
