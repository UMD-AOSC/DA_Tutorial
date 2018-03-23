from class_maooam import maooam
import numpy as np
from class_state_vector import state_vector

import params_maooam
from params_maooam import ndim, tw, t_run, t_trans, dt
#from maooam import params_maooam
#from maooam.params_maooam import ndim, tw, t_run, t_trans, dt
from maooam import integrator
import time
import sys

#------------------------------------------------------------------
# Setup initial state
#------------------------------------------------------------------
ic_file='x0.dat'
name = 'x_nature'

# Get definition of time from params_maooam.py file:
t0=0.0
tf=t_run
dt=dt
if t_trans > 0:
  trans_tvec = np.arange(t0-t_trans/dt, t_trans, dt)
tvec = np.arange(t0, tf, dt)
state0=np.zeros(ndim)

#------------------------------------------------------------------
# (Optional) Update the initial starting point
#------------------------------------------------------------------
 # Load initial state
try:
  with open(ic_file) as f:
    lines = f.readlines()
  x0a = np.array(lines)
  state0 = x0a.astype(np.float)
except:
  print('File does not exist: ', ic_file)
  print('To specify initial conditions, add file: ', ic_file)

print ('x0 = ', state0)
params = []

#------------------------------------------------------------------
# (Optional) Add initial perturbation
#------------------------------------------------------------------
if len(sys.argv) > 1 and sys.argv[1] == "freerun":
  name = 'x_freerun'
  initial_perturbation = np.squeeze(0.01*(np.random.rand(1,ndim)*2-1))
  print('initial_perturbation = ', initial_perturbation)
  state0 = state0 + initial_perturbation
  print('initial state = ', state0)

#------------------------------------------------------------------
# Setup state vector object
#------------------------------------------------------------------
sv = state_vector(params=params,x0=state0,t=tvec,name=name)

#------------------------------------------------------------------
# Initialize the maooam object
#------------------------------------------------------------------
model = maooam()

#------------------------------------------------------------------
# Run MAOOAM to generate a nature run with the specified parameters
#------------------------------------------------------------------
print('Run MAOOAM...')
if t_trans > 0:
  trans_trajectory = model.run(sv.x0,sv.t)
  sv.x0 = trans_trajectory[-1,:]
  
trajectory = model.run(sv.x0,sv.t)
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

