#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

import pandas
import sys
import os
DIR = os.path.join(os.path.dirname(__file__), "ref")

field_name = sys.argv[1]
ref = sys.argv[2]

# from pipe_params import *

A = pandas.read_csv("plot.csv")

v0 = 0.5

rs = A["Points:0"]
vys = A[field_name + ":1"]

# ref data
if ref == "vdi-ddp10-re1":
    A = pandas.read_csv(os.path.join(DIR, "velocity-profile-VDI-M7-Re-1.csv"))
    ddp = 10
elif ref == "witso-ddp4-re1":
    A = pandas.read_csv(os.path.join(DIR, "WiTso2000-eta_f-Re1.csv"))
    ddp = 4
elif ref == "witso-ddp20-re1":
    A = pandas.read_csv(os.path.join(DIR, "WiTso2000-Re1-Ddp-20.csv"))
    ddp = 20
elif ref == "witso-ddp4-eff-re1":
    A = pandas.read_csv(os.path.join(DIR, "WiTso2000-Re1-eta_eff.csv"))
    ddp = 4
elif ref == "witso-ddp4-re1000":
    A = pandas.read_csv(os.path.join(DIR, "WiTso2000-eta_f-Re1000.csv"))
    ddp = 4
else:
    assert False

bed_radius = max(rs)
print("bed radius is", bed_radius)
pellet_diameter = 2*bed_radius / ddp

vys_scaled = -vys / v0
rs_scaled = (bed_radius - rs)/pellet_diameter

rs_ref = A["x"]
vys_ref = A["v"]

fig, ax = plt.subplots()

ax.plot(rs_scaled, vys_scaled, label="num")
ax.plot(rs_ref, vys_ref, label="ref")

ax.set_xlim(0,2)
# ax.set_ylim(0,4)
ax.axhline(1)
ax.legend()

fig.savefig("plot.png")

plt.close(fig)
