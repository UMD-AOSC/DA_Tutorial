#!/bin/python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors, ticker, cm
from class_state_vector import state_vector

name = 'x_nature'
infile = name+'_LEs.pkl'
sv = state_vector()
sv = sv.load(infile)
#LEs_trans = np.load('LEs_transient.dat.npy')

lyap_exp = sv.getLEs()

print (lyap_exp)

t1 = np.arange(-1.0, 7.0, 1.0)

plt.plot(t1,np.zeros_like(t1),'k--',linewidth=1.0)
plt.plot(lyap_exp, 'ro')

plt.xlim([-0.5,3.5])

plt.show()
