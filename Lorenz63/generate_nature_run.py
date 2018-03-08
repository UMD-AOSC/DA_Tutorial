from class_lorenz63 import lorenz63
import numpy as np
from class_state_vector import state_vector

#------------------------------------------------------------------
# Setup initial state
#------------------------------------------------------------------
sigma = 10.0
rho = 28.0
beta = 8.0 / 3.0
params = [sigma, rho, beta]

name = 'x_nature'

t0=0.0
tf=40.0
dt=0.001
state0 = [1.0, 1.0, 1.0]
tvec = np.arange(t0, tf, dt)

#------------------------------------------------------------------
# (Optional) Update the initial starting point
#------------------------------------------------------------------
#state0 = [-1.6831903,  -2.78556131, 13.0107312]
#state0 = [100,100,100]

#------------------------------------------------------------------
# (Optional) Add initial perturbation
#------------------------------------------------------------------
#name = 'x_freerun'
#initial_perturbation = np.squeeze(0.01*(np.random.rand(1,3)*2-1))
#print('initial_perturbation = ', initial_perturbation)
#climate_std =  [7.44085386, 8.38759537, 7.84386348]
#print('climate_std = ', climate_std)
#state0 = state0 + initial_perturbation*climate_std
#print('initial state = ', state0)

#------------------------------------------------------------------
# Setup state vector object
#------------------------------------------------------------------
sv = state_vector(params=params,x0=state0,t=tvec,name=name)

#------------------------------------------------------------------
# Initialize the l63 object
#------------------------------------------------------------------
l63 = lorenz63(sigma=sigma, rho=rho, beta=beta)

#------------------------------------------------------------------
# Run L63 to generate a nature run with the specified parameters
#------------------------------------------------------------------
print('Run l63...')
trajectory = l63.run(sv.x0,sv.t)
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

