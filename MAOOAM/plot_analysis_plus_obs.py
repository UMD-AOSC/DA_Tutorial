from class_maooam import maooam
from class_state_vector import state_vector
from class_obs_data import obs_data
from class_da_system import da_system
from sys import argv

#------------------------------------------------------------------
# Read state vector object
#------------------------------------------------------------------
name = 'x_analysis'

#----------------
# Choose method:
#----------------
method = argv[1]
infile=name+'_'+method+'.pkl'
das = da_system()
das = das.load(infile)

sv = das.getStateVector()
obs = das.getObsData()

print(das)
print(sv)
print(obs)

#------------------------------------------------------------------
# Plot the result
#------------------------------------------------------------------
model = maooam()
title = 'MAOOAM attractor using DA method: %s'%(method)

states = sv.getTrajectory()
points=obs.getVal()
error = abs(obs.getVal()-sv.getTrajectory())

print('states = ')
print(states)

print('points = ')
print(points)

print('difference = ')
print(error)

model.plot_lines_and_points(states=states,points=points,cvec=error,plot_title=title)
