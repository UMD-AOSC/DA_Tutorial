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
infile = name+'_Mhist.pkl'
sv = state_vector()
sv = sv.load(infile)
x_nature = sv.getTrajectory()
t_nature = sv.getTimes()

#-----------------------------------------------------------------------
# Initialize the timesteps
#-----------------------------------------------------------------------
t_nature = sv.getTimes()
ainc_step = 10  # (how frequently to perform an analysis)
dtau = (t_nature[ainc_step] - t_nature[0])
maxit,xdim = np.shape(x_nature)

#------------------------------------------------------------------
# Compute the LEs by propagating the TLM in time and intermittently 
# performing a QR decomposition and rescaling as necessary.
#------------------------------------------------------------------
I = np.identity(xdim)
Mhist = sv.getMhist()
#print('Mhist = ')
#print(Mhist)
rescale_interval = ainc_step

M2hist = []
Qhist = []
Rhist = []
hist_idx = []
M2 = I
for i in range(maxit):
  
  # Iterate the linear propagator (evolution matrix) 
  # to fit the specified rescaling time interval

# print('M = ')
# print(Mhist[i])
  M2 = np.dot(Mhist[i],M2)

  # Rescale via QR decomposition
  if (np.mod(i+1,rescale_interval)==0):
#   print('=======================================================')
#   print('M2 = ')
#   print(M2)

    Q,R = np.linalg.qr(M2)

#   print('Q = ')
#   print(Q)
#   print('R = ')
#   print(R)
    
    # Store the Q and R matrices
    M2hist.append(deepcopy(M2))
    Qhist.append(deepcopy(Q))
    Rhist.append(deepcopy(R))
    hist_idx.append(i)
    M2 = Q

print('Last Q = ')
print(Q)

print('Last R = ')
print(R)

sv.rescale_interval = rescale_interval
sv.hist_ainc = ainc_step
sv.hist_dtau = dtau
sv.hist_idx = hist_idx
sv.setM2hist(M2hist)
sv.setQhist(Qhist)
sv.setRhist(Rhist)

#------------------------------------------------------------------
# Store the nature run data
#------------------------------------------------------------------
outfile = name+'_QR.pkl'
sv.save(outfile)
