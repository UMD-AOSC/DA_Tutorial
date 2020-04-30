from class_lorenz96 import lorenz96
import numpy as np
from class_state_vector import state_vector

#------------------------------------------------------------------
# Setup initial state
#------------------------------------------------------------------
xdim = 5
F = 8.0
params = [F]

name = 'x_nature'

t0=0.0
tf=100.0
dt=0.001
state0 = np.ones(xdim) + np.random.uniform(low=-1.0, high=1.0, size = (xdim))
tvec = np.arange(t0, tf, dt)

#------------------------------------------------------------------
# (Optional) Update the initial starting point
#------------------------------------------------------------------
#state0 = [-1.6831903,  -2.78556131, 13.0107312, 1, 1]
#state0 = [100,100,100,100,100]

#------------------------------------------------------------------
# (Optional) Add initial perturbation
#------------------------------------------------------------------
# From previous run:
#Climatological Mean:
#[2.37742045 2.26356092 2.45364065 2.31241994 2.02554909]
#Climatological Standard Deviation:
#[3.63477096 3.63919927 3.30788215 3.76514026 3.62503822]

#name = 'x_freerun'
#initial_perturbation = np.squeeze(0.01*(np.random.rand(1,3)*2-1))
#print('initial_perturbation = ', initial_perturbation)
#climate_std =  [3.63477096 3.63919927 3.30788215 3.76514026 3.62503822]
#print('climate_std = ', climate_std)
#state0 = state0 + initial_perturbation*climate_std
#print('initial state = ', state0)

#------------------------------------------------------------------
# Setup state vector object
#------------------------------------------------------------------
sv = state_vector(params=params,x0=state0,t=tvec,name=name)

#------------------------------------------------------------------
# Initialize the l96 object
#------------------------------------------------------------------
l96 = lorenz96(F=F)

#------------------------------------------------------------------
# Run L96 to generate a nature run with the specified parameters
#------------------------------------------------------------------
print('Run l96...')
trajectory = l96.run(sv.x0,sv.t)
sv.setTrajectory(trajectory)

#------------------------------------------------------------------
# Output the beginning and end states, and compute model climatology
#------------------------------------------------------------------
print(sv)

#------------------------------------------------------------------
# Store the nature run data
#------------------------------------------------------------------
outfile = name+'.pkl'
sv.save(outfile)

