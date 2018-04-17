#!/usr/bin/env python
from ctypes import *
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

n = 36

class MaooamFortran:
    module_maooam = np.ctypeslib.load_library("step_maooam.so", ".")
    module_maooam.step_maooam_.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64)]
    module_maooam.step_maooam_.restype = c_void_p

    def __init__(self, dt):
        assert dt.__class__ in [float, np.float32, np.float64]
        self.dt = np.array([dt])

    def step(self, x0):
        self.module_maooam.step_maooam_(x0, self.dt)
        return x0

def main():
    nt = 10000000
    intvl = 100
    mf = MaooamFortran(0.01)
    x0 = __model_state_example()
    xhist = np.empty((nt // intvl, n))
    for i in range(nt):
        x0 = mf.step(x0)
        if i % intvl == 0:
            xhist[i // intvl, :] = x0
    print(x0)
    plot_xhist(xhist, "python.pdf")
    plot_fortran()

def plot_xhist(data, imgname):
    cm = plt.imshow(data, aspect="auto")
    cm.set_clim(-0.03, 0.03)
    plt.colorbar(cm)
    plt.savefig(imgname)
    plt.close()

def __model_state_example():
    xini = np.array([
        4.695340259215241E-002,
        2.795833230987369E-002,
        -2.471191763590483E-002,
        -7.877635082773315E-003,
        -4.448292568544942E-003,
        -2.756238610924190E-002,
        -4.224400051368891E-003,
        5.914241112882518E-003,
        -1.779437742222920E-004,
        5.224450720394076E-003,
        4.697982667229096E-002,
        5.149282577209392E-003,
        -1.949084549066326E-002,
        4.224006062949761E-004,
        -1.247786759371923E-002,
        -9.825952138046594E-003,
        -2.610941795170075E-005,
        2.239286581216401E-003,
        -7.891896725509534E-004,
        7.470171905055880E-004,
        -9.315932162526787E-007,
        3.650179005106874E-005,
        1.064122403269511E-006,
        3.937836448211443E-008,
        -2.208288760403859E-007,
        -3.753762121228048E-006,
        -7.105126469908465E-006,
        1.518110190916469E-008,
        -5.773178576933025E-004,
        0.187369278208256,
        1.369868543156558E-003,
        7.023608700166264E-002,
        -4.539810680860224E-004,
        -1.882650440363933E-003,
        -3.900412687995408E-005,
        -1.753655087903711E-007])
    return xini

def read_text_dat(filename):
    with open(filename, "r") as f:
        # lines = f.readlines()
        words = f.read().split()
    lines = [words[i:i + n + 1] for i in range(0, len(words), n + 1)]
    nl = len(lines)
    xhist = np.empty((nl, n))
    for i, l in enumerate(lines):
        li = list(map(float, l))
        xhist[i, :] = li[1:]
    return xhist

def plot_fortran():
    xhist = read_text_dat("evol_field.dat")
    plot_xhist(xhist, "fortran.pdf")

if __name__ == "__main__":
    main()

