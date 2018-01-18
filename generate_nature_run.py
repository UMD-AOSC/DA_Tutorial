from class_lorenz63 import lorenz63
import numpy as np
from class_state_vector import state_vector

#------------------------------------------------------------------
# Setup initial state
#------------------------------------------------------------------
rho = 28.0
sigma = 10.0
beta = 8.0 / 3.0
params = [sigma, rho, beta]

outfile = 'x_nature.pkl'

t0=0.0
tf=40.0
dt=0.01
state0 = [1.0, 1.0, 1.0]
tvec = np.arange(t0, tf, dt)

#------------------------------------------------------------------
# Setup state vector object
#------------------------------------------------------------------
sv = state_vector(params=params,x0=state0,t=tvec,name='x_nature')

#------------------------------------------------------------------
# Initialize the l63 object
#------------------------------------------------------------------
l63 = lorenz63(sigma=sigma, rho=rho, beta=beta)

#------------------------------------------------------------------
# Run L63 to generate a nature run
#------------------------------------------------------------------
print('Run l63...')
sv.trajectory = l63.run(sv.x0,sv.t)

#------------------------------------------------------------------
# Store the nature run data
#------------------------------------------------------------------
sv.save(outfile)

