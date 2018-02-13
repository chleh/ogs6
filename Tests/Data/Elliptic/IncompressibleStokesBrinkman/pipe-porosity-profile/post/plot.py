#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

import pandas

A = pandas.read_csv("plot.csv")

v0 = 0.5
r_bed = 1
d_pel = 0.5

rs = A["Points:0"]
vys = A["darcy_velocity:1"]

vys_scaled = -vys / v0
rs_scaled = (r_bed - rs)/d_pel


# ref data
if False:
    A = pandas.read_csv("velocity-profile-VDI-M7-Re-1.csv")
elif False:
    A = pandas.read_csv("WiTso2000-eta_f-Re1.csv")
else:
    A = pandas.read_csv("WiTso2000-eta_f-Re1000.csv")
rs_ref = A["x"]
vys_ref = A["v"]

fig, ax = plt.subplots()

ax.plot(rs_scaled, vys_scaled)
ax.plot(rs_ref, vys_ref)

ax.set_xlim(0,2)
ax.set_ylim(0,4)
ax.axhline(1)

fig.savefig("plot.png")

plt.close(fig)


fig, ax = plt.subplots()

drs = np.diff(rs_scaled)
dvys = np.abs(np.diff(vys_scaled))

ax.plot(rs_scaled[1:], dvys/drs)

fig.savefig("plot-grad-vy.png")

plt.close(fig)

if rs_scaled[0] != r_bed/d_pel:
    print(rs_scaled[0], r_bed/d_pel)
    print(rs_scaled)

print(np.trapz(x=rs_scaled, y=(r_bed - rs_scaled*d_pel)  * vys_scaled))
