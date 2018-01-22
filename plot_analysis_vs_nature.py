from class_lorenz63 import lorenz63
from class_state_vector import state_vector
from class_obs_data import obs_data
from class_da_system import da_system
import numpy as np

#------------------------------------------------------------------
# Read state vector objects
#------------------------------------------------------------------
name = 'x_analysis'
method = 'ETKF'
infile1 = name+'_'+method+'.pkl'
das = da_system()
das = das.load(infile1)
sv1 = das.getStateVector()
print(sv1)

infile2='x_nature.pkl'
sv2 = state_vector()
sv2 = sv2.load(infile2)
print(sv2)

#------------------------------------------------------------------
# Plot the result
#------------------------------------------------------------------
error = np.abs(sv2.getTrajectory() - sv1.getTrajectory())
times = sv1.getTimes()
l63 = lorenz63()
title='Comparison of '+method+' analysis and nature run'
l63.plot_lines_and_lines(states1=sv1.getTrajectory(),states2=sv2.getTrajectory(),cvec=times,name1=sv1.name,name2=sv2.name,plot_title=title)
