"""
    Principal module
    ======================

    Python implementation of the Modular Arbitrary-Order Ocean-Atmosphere Model MAOOAM

    .. note :: The python code is available here : \
    `maooam.py <../_modules/maooam.html>`_ .


    :Example:

    >>> from maooam import *

    Global variable
    -------------------

    * **ic.X0** : initial conditions
    * **X** : live step vector
    * **t** : time
    * **t_trans**, **t_run** : respectively transient and running time
    * **dt** : step time
    * **tw** : step time for writing on evol_field.dat

    Dependencies
    -------------------

    >>> import numpy as np
    >>> import params_maooam
    >>> from params_maooam import ndim,tw,t_run,t_trans,dt
    >>> import aotensor
    >>> import time
    >>> import ic_def
    >>> import ic
    >>> import sys
"""

import numpy as np
import params_maooam
from params_maooam import ndim, tw, t_run, t_trans, dt
import integrator
import time
import ic_def
import ic
import sys

def print_progress(p):
    sys.stdout.write('Progress {:.2%} \r'.format(p))
    sys.stdout.flush()

class bcolors:
    """to color the instructions in the console"""
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


print (bcolors.OKBLUE + "Model MAOOAM v1.3" + bcolors.ENDC)
print (bcolors.OKBLUE + "Initialization ..." + bcolors.ENDC)

ic_def.load_IC()

X = ic.X0
print (bcolors.OKBLUE + "Starting the transient time evolution ..." + bcolors.ENDC)
t = 0.
T = time.clock()
t_up = dt/t_trans*100
while t < t_trans:
    X = integrator.step(X, t, dt)
    t += dt
    if t/t_trans*100 % 0.1 < t_up:
        print_progress(t/t_trans)

print (bcolors.OKBLUE + "Starting the time evolution ..." + bcolors.ENDC)
fichier = open("evol_field.dat", "w")
t = 0.
t_up = dt/t_run*100

while t < t_run:
    X = integrator.step(X, t, dt)
    t += dt
    if t % (tw) < dt:
        fichier.write(str(t)+" ")
        for i in range(0, ndim):
            fichier.write(str(X[i])+" ")
        fichier.write("\n")
    if t/t_run*100 % 0.1 < t_up:
        print_progress(t/t_run)
fichier.close()
print (bcolors.OKBLUE + "Evolution finished " + bcolors.ENDC)

print (bcolors.OKBLUE + "Time clock :" + bcolors.ENDC)
print (time.clock()-T)

