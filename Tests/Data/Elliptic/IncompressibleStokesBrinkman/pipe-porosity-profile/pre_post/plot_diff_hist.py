#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

for fn in \
    "diff_hist_abs_diff_comp_0.npz", \
    "diff_hist_abs_diff_comp_1.npz", \
    "diff_hist_abs_diff_comp_2.npz", \
    "diff_hist_rel_diff_comp_0.npz", \
    "diff_hist_rel_diff_comp_1.npz", \
    "diff_hist_rel_diff_comp_2.npz":

    try:
        A = np.load(fn)
    except FileNotFoundError:
        continue

    hist = A["hist"]
    edges = A["edges"]

    fig, (ax, ax2) = plt.subplots(2, sharex=True)

    widths = np.diff(edges)
    centers = edges[:-1] + 0.5*widths

    # ax.bar(centers, hist, 0.8*widths)
    ax.plot(centers, hist, color="r")
    ax.set_xlim(edges[0], edges[-1])
    ax.set_ylim(0.5)
    ax.set_yscale("log")

    ax2.plot(centers, hist, color="r")

    fig.savefig(fn[:-len("npz")] + "png")
