#!/usr/bin/env python3

import sys
import subprocess as sp
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NDIM = 36

def plot_time_colormap(dat, img_name, vmin=None, vmax=None, title="", cmap="RdBu_r"):
    assert dat.__class__ == np.ndarray
    assert len(dat.shape) == 2
    assert dat.shape[1] == NDIM
    cm = plt.imshow(dat, aspect="auto", cmap=cmap)
    if (vmin is not None) and (vmax is not None):
        cm.set_clim(vmin, vmax)
    plt.colorbar(cm)
    plt.xlabel("model variable")
    plt.ylabel("time")
    plt.title(title)
    plt.savefig(img_name)
    plt.close()

def __test_plot_time_colormap():
    nt = 100
    dat = np.random.randn(nt, NDIM)
    sp.run("mkdir -p img", shell=True, check=True)
    plot_time_colormap(dat, "img/tmp.png", None, None, "test")

if __name__ == "__main__":
    __test_plot_time_colormap()
