from class_lorenz63 import lorenz63
from class_state_vector import state_vector
from class_obs_data import obs_data
from class_da_system import da_system
import numpy as np
from sys import argv

#------------------------------------------------------------------
# Read state vector objects
#------------------------------------------------------------------
name = 'x_analysis'

method1 = argv[1]
infile1 = name+'_'+method1+'.pkl'
das = da_system()
das = das.load(infile1)
sv1 = das.getStateVector()
print(sv1)

method2 = argv[2]
infile2 = name+'_'+method2+'.pkl'
das = da_system()
das = das.load(infile2)
sv2 = das.getStateVector()
print(sv2)

#------------------------------------------------------------------
# Plot the result
#------------------------------------------------------------------
outfile=method1+'_vs_'+method2
error = np.abs(sv2.getTrajectory() - sv1.getTrajectory())
times = sv1.getTimes()
l63 = lorenz63()
title='Comparison of '+method1+' analysis and '+method2+' analysis'
l63.plot_lines_and_lines(states1=sv1.getTrajectory(),states2=sv2.getTrajectory(),cvec=times,name1=sv1.name,name2=sv2.name,plot_title=title,outfile=outfile)
