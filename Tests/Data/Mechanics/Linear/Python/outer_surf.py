#!/usr/bin/python

import numpy as np

CELLS_PER_DIMENSION = 100
N = 2*CELLS_PER_DIMENSION + 1


for i in range(N):
    phi = np.pi * 0.5 / (N-1) * i
    x = np.sin(phi)
    y = np.cos(phi)
    print(f'<point id="{i+3}" x="{x}" y="{y}" z="0" />')

for i in range(N):
    print(f'<pnt>{i+3}</pnt>')
