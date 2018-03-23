#!/usr/bin/env python3

import sys
import subprocess as sp
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from class_state_vector import state_vector
from class_da_system import da_system

NDIM = 36

def plot_time_colormap(dat, img_name, vmin=None, vmax=None, title="", cmap="RdBu_r"):
    assert dat.__class__ == np.ndarray
    assert len(dat.shape) == 2
    assert dat.shape[1] == NDIM
    cm = plt.imshow(dat, aspect="auto", cmap=cmap, origin="bottom")
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

def __sample_read_files():
    vlim_raw = [-0.05, 0.1]
    vlim_diff = [None, None]
    nature_file ='x_nature.pkl'
    nature = state_vector()
    nature = nature.load(nature_file)
    freerun_file = 'x_freerun.pkl'
    freerun = state_vector()
    freerun = freerun.load(freerun_file)
    sp.run("mkdir -p img", shell=True, check=True)
    plot_time_colormap(freerun.getTrajectory() - nature.getTrajectory(),
                       "img/error_free_run.png", *vlim_diff, "error free run")
    plot_time_colormap(freerun.getTrajectory(),
                       "img/freerun.png", *vlim_raw, "freerun", "viridis")
    plot_time_colormap(nature.getTrajectory(),
                       "img/nature.png", *vlim_raw, "nature", "viridis")
    for method in ["ETKF", "nudging"]:
        analysis_file = 'x_analysis_{method}.pkl'.format(method=method)
        das = da_system()
        das = das.load(analysis_file)
        analysis = das.getStateVector()
        plot_time_colormap(analysis.getTrajectory() - nature.getTrajectory(),
                           "img/error_analysis_%s.png" % method, *vlim_diff,
                           "error analysis %s" % method)
        plot_time_colormap(analysis.getTrajectory(),
                           "img/analysis_%s.png" % method, *vlim_raw,
                           "analysis %s" % method, "viridis")

if __name__ == "__main__":
    __sample_read_files()
