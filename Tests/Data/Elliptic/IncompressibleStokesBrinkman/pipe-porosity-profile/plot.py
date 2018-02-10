#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

import pandas

# A = np.loadtxt("plot.csv", unpack=True)
A = pandas.read_csv("plot.csv") #, unpack=True)

v0 = 0.5
r_bed = 1
d_pel = 0.2

rs = A["Points:0"]
vys = A["darcy_velocity:1"]

vys_scaled = -vys / v0
rs_scaled = (r_bed - rs)/d_pel


# ref data
A = pandas.read_csv("velocity-profile-VDI-M7-Re-1.csv")
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
