#!/usr/bin/env pvpython

from __future__ import print_function

import paraview.simple as pview
import cl_work.pre.via_paraview as pre

import numpy as np
# import scipy as sp
from scipy.integrate import odeint

cp = 4000

MV = 0.01806
MN = 0.028013

phi = 0.5

def MG(xmV):
    return MV * MN / (MN*xmV + MV*(1.0-xmV))

def dMGdx_over_MG(xmV):
    return (MV - MN) / (MN*xmV + MV*(1.0-xmV))

p0 = 1e5
T0 = 3e2
xmV0 = 0.5
R = 8.3144598

rhoSR0 = 1
rhoGR0 = p0 * MG(xmV0) / R / T0
rho0 = phi * rhoGR0 + (1-phi)*rhoSR0

rho_hat = (1-phi) * 0.1
h_ads = 4e6

def f(x, t):
    p, T, xmV = x
    rhoG = phi * p * MG(xmV) / R / T

    A = np.matrix([
        -phi, cp * rho0, 0,
        rhoG/p, -rhoG/T, rhoG * dMGdx_over_MG(xmV),
        0, 0, rhoG
        ]).reshape((3,3))

    b = np.array([ rho_hat*h_ads, -rho_hat, -(1.0-xmV)*rho_hat ]).reshape((3,1))
    return np.linalg.solve(A, b).reshape(3)


print(f([p0, T0, xmV0], 0))

ts = np.linspace(0, 1, 101)
res = odeint(f, [p0, T0, xmV0], ts, rtol=1e-12, atol=1e-12)

ps = res[:,0]
Ts = res[:,1]
xs = res[:,2]

if False:
    import matplotlib.pyplot as plt
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True)

    ax1.plot(ts, ps)
    ax2.plot(ts, Ts)
    ax3.plot(ts, xs)

    rhoS = rhoSR0 * phi + ts * rho_hat
    rhoG = (1.0 - phi) * ps * MG(xs) / R / Ts
    ax4.plot(ts, rhoS)
    ax4.plot(ts, rhoG)
    ax4.plot(ts, rhoGR0 * (1-phi) - ts * rho_hat)
    ax4.plot(ts, rhoS + rhoG)

    fig.savefig("reference_results.png")


i = 0

def f_profile_T(coords):
    import numpy as np
    return Ts[i] * np.ones(coords.shape[0])

def f_p_ref(coords):
    import numpy as np
    return ps[i] * np.ones(coords.shape[0])

def f_xmV_ref(coords):
    import numpy as np
    return xs[i] * np.ones(coords.shape[0])

def f_v_Darcy_ref(coords):
    import numpy as np
    return np.zeros((coords.shape[0], 2))

def f_rate_ref(coords):
    import numpy as np
    return np.ones(coords.shape[0]) * 0.1

def f_density_ref(coords):
    import numpy as np
    return np.ones(coords.shape[0]) * (1.0 + ts[i] * 0.1)

reader = pview.XMLUnstructuredGridReader(FileName=["pipe_quadratic.vtu"])

for i in np.where(10*ts % 2 == 0)[0]:
    print(i, ts[i], ps[i], Ts[i], xs[i])
    T_ref = pre.NodalFunction("T_ref", f_profile_T, Input=reader)
    p_ref = pre.NodalFunction("p_ref", f_p_ref, Input=T_ref)
    p_ref.CopyArrays = 1
    xmV_ref = pre.NodalFunction("xmV_ref", f_xmV_ref, Input=p_ref)
    xmV_ref.CopyArrays = 1
    v_ref = pre.NodalFunction("v_Darcy_ref", f_v_Darcy_ref, Input=xmV_ref)
    v_ref.CopyArrays = 1
    rate_ref = pre.NodalFunction("reaction_rate_ref", f_rate_ref, Input=v_ref)
    rate_ref.CopyArrays = 1
    dens_ref = pre.NodalFunction("solid_density_ref", f_density_ref, Input=rate_ref)
    dens_ref.CopyArrays = 1

    # save data
    pview.SaveData("pipe_ref_{}s.vtu".format(ts[i]), proxy=dens_ref, DataMode='Binary',
        EncodeAppendedData=1,
        CompressorType='ZLib')
