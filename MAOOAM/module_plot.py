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

def __sample_read_files():
    for method in ["nudging", "ETKF"]:
        nature_file ='x_nature.pkl'
        nature = state_vector()
        nature = nature.load(nature_file)
        freerun_file = 'x_freerun.pkl'
        freerun = state_vector()
        freerun = freerun.load(freerun_file)
        analysis_file = 'x_analysis_{method}.pkl'.format(method=method)
        das = da_system()
        das = das.load(analysis_file)
        analysis = das.getStateVector()
        plt.plot(nature.getTimes(),
                 np.linalg.norm(freerun.getTrajectory() - nature.getTrajectory(),
                                axis=1), label='Free run')
        plt.plot(nature.getTimes(),
                 np.linalg.norm(analysis.getTrajectory() - nature.getTrajectory(),
                                axis=1), label='Analysis ({method})'.format(method=method))
        plt.legend()
        plt.xlabel('Time')
        plt.ylabel('Error', rotation='horizontal', labelpad=20)
        plt.tight_layout()
        plt.savefig("img/tmp.png")
        plt.close()

if __name__ == "__main__":
    __sample_read_files()
