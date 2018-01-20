# Tutorial: "A tour of Data Assimilation methods"
# Model: Lorenz-63
# DA Methods: Nudging, 3D-Var, 4D-Var, Particle Filter, EnKF, Hybrid
import numpy as np
from class_lorenz63 import lorenz63
from class_state_vector import state_vector
from class_obs_data import obs_data
from class_da_system import da_system
from copy import deepcopy

#-----------------------------------------------------------------------
# Read the L63 nature run
#-----------------------------------------------------------------------
name = 'x_nature'
infile = name+'_Jhist.pkl'
sv = state_vector()
sv = sv.load(infile)
x_nature = sv.getTrajectory()
t_nature = sv.getTimes()

#-----------------------------------------------------------------------
# Initialize the timesteps
#-----------------------------------------------------------------------
t_nature = sv.getTimes()
ainc_step = 4  # (how frequently to perform an analysis)
dtau = (t_nature[ainc_step] - t_nature[0])
tsteps=10 * ainc_step
dt = dtau/tsteps
maxit,xdim = np.shape(x_nature)

#------------------------------------------------------------------
# Compute the LEs by propagating the TLM in time and intermittently 
# performing a QR decomposition and rescaling as necessary.
#------------------------------------------------------------------
I = np.identity(xdim)
Jhist = sv.getJhist()
print('Jhist = ')
print(Jhist)
rescale_interval = ainc_step

Mhist = []
Qhist = []
Rhist = []
hist_idx = []
M = I
for i in range(maxit-1):
  
  # Iterate the linear propagator (evolution matrix) 
  # to fit the specified rescaling time interval

# print('J = ')
# print(Jhist[i])
  M = np.dot(Jhist[i],M)

  # Rescale via QR decomposition
  if (np.mod(i+1,rescale_interval)==0):
    print('=======================================================')
    print('M = ')
    print(M)

    Q,R = np.linalg.qr(M)

    print('Q = ')
    print(Q)
    print('R = ')
    print(R)
    
    # Store the Q and R matrices
    Mhist.append(deepcopy(M))
    Qhist.append(deepcopy(Q))
    Rhist.append(deepcopy(R))
    hist_idx.append(i)
    M = Q

sv.rescale_interval = rescale_interval
sv.hist_idx = hist_idx
sv.setMhist(Mhist)
sv.setQhist(Qhist)
sv.setRhist(Rhist)

#------------------------------------------------------------------
# Store the nature run data
#------------------------------------------------------------------
outfile = name+'_TLM.pkl'
sv.save(outfile)
