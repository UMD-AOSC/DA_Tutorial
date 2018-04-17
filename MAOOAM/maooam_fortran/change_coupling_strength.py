#!/usr/bin/env python3

import sys
import subprocess as sp
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def update_params_nml(coupl):
    assert coupl.__class__ in [int, float, np.float32, np.float64]
    fname = "params.nml"
    line_num = 19
    with open(fname) as f:
        st = f.read().splitlines()
    st2 = []
    for i, l in enumerate(st):
        if i + 1 == line_num:
            assert "COUPLE = " in l
            l = "  COUPLE = %f" % coupl
        st2.append(l)
    st2 = "\n".join(st2)
    with open(fname, "w") as f:
        f.write(st2)

def update_plot_tlm_py(coupl):
    assert coupl.__class__ in [int, float, np.float32, np.float64]
    fname = "plot_tlm.py"
    line_num = 16
    with open(fname) as f:
        st = f.read().splitlines()
    st2 = []
    for i, l in enumerate(st):
        if i + 1 == line_num:
            assert "FNAME = " in l
            l = "FNAME = \"dat/evol_field_tlm_%.2f.dat\"" % coupl
        st2.append(l)
    st2 = "\n".join(st2)
    with open(fname, "w") as f:
        f.write(st2)

def plot_le(coupls):
    n = len(coupls)
    for i in range(n):
        fname = "dat/le_%.2f.txt" % coupls[i]
        with open(fname) as f:
            st = f.read().splitlines()
        l = st[1]
        l = l.replace("[", "").replace("]", "")
        les = np.array(list(map(float, l.split())))
        ones = np.ones(len(les))
        plt.scatter(ones * coupls[i], les, marker=".", color="black", alpha=0.3)
    sp.run("mkdir -p img", shell=True, check=True)
    plt.axhline(0.0, zorder=-1, alpha=0.5)
    plt.xlabel("coupling strength relative to VL2016")
    plt.ylabel("lyapunov exponents")
    plt.savefig("img/tmp.png")

def loop():
    n = 21
    coupls = np.linspace(0.0, 2.0, n)
    sp.run("make")
    sp.run("mkdir -p dat", shell=True, check=True)
    for i in range(n):
        update_params_nml(coupls[i])
        sp.run("./run_true_and_tlm", shell=True, check=True)
        sp.run("mv evol_field_tlm.dat dat/evol_field_tlm_%.2f.dat" % coupls[i], shell=True, check=True)
        update_plot_tlm_py(coupls[i])
        sp.run("python3 plot_tlm.py > dat/le_%.2f.txt" % coupls[i], shell=True, check=True)
    plot_le(coupls)
    
if __name__ == "__main__":
    loop()
    
