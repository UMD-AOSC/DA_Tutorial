# Tutorial: "A tour of Data Assimilation methods"
# Model: Lorenz-63
# DA Methods: Nudging, 3D-Var, 4D-Var, Particle Filter, EnKF, Hybrid
import numpy as np
from class_lorenz63 import lorenz63
from class_state_vector import state_vector
from class_obs_data import obs_data
from class_da_system import da_system
from sys import argv
from copy import deepcopy

#-----------------------------------------------------------------------
# Usage:
#  python generate_analysis_LEs.py <method>
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
# Read the L63 nature run
#-----------------------------------------------------------------------
name1 = 'x_nature'
infile1 = name1+'_Mhist.pkl'
sv1 = state_vector()
sv1 = sv1.load(infile1)
x_nature = sv1.getTrajectory()
t_nature = sv1.getTimes()

#-----------------------------------------------------------------------
# Read the analysis
#-----------------------------------------------------------------------
name = 'x_analysis'
method = argv[1]
infile2 = name+'_'+method+'.pkl'
das = da_system()
das = das.load(infile2)
sv = das.getStateVector()
KH, khidx = das.getKH()
print(sv)

#print('KH = ')
#print(KH)
#print('khidx = ')
#print(khidx)

#-----------------------------------------------------------------------
# Initialize the timesteps
#-----------------------------------------------------------------------
t_nature = sv1.getTimes()
dt = (t_nature[1] - t_nature[0])
_,xdim = np.shape(x_nature)
acyc_step = das.acyc_step
dtau = das.dtau

#------------------------------------------------------------------
# Propagate the TLM in time and intermittently 
# performing a QR decomposition to rescale as necessary.
#------------------------------------------------------------------
rescale_interval = 2*acyc_step
I = np.identity(xdim)
Mhist = sv1.getMhist()
M2hist = []
Qhist = []
Rhist = []
hist_idx = []
M2 = I
maxit = das.maxit
ki=0
for i in range(maxit - acyc_step):
  
  # In order to find the eigenvalues of the final matrix describing
  # the evolution of the entire trajectory, we would like to 
  # compute the product of the M matrices to identify the
  # stretching and contracting that results.
  #
  # Iterate the linear propagator (evolution matrix) 
  # to fit the specified rescaling time interval

  M2 = np.dot(Mhist[i],M2)

  # If there was an analysis at this time index, then
  # apply the contraction (I-KH) applied by the DA
  if khidx[ki] == i:
    M2 = np.dot(I-KH[ki],M2) 
    ki = ki+1

  # Rescale via QR decomposition
  if (np.mod(i+1,rescale_interval)==0):
    Q,R = np.linalg.qr(M2)
    
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
sv.hist_ainc = acyc_step
sv.hist_dtau = dtau
sv.hist_idx = hist_idx
sv.setM2hist(M2hist)
sv.setQhist(Qhist)
sv.setRhist(Rhist)

#-----------------------------------------------------------------------
# Compute the Lyapunov Exponents (LEs) of the data assimilation system
#-----------------------------------------------------------------------
lyap_sum = np.zeros(xdim)
rescale_cnt = 0
for i in range(len(Rhist)):
  R = Rhist[i]
  rdiag = abs(np.diag(R))
  lyap_sum = lyap_sum + np.log(rdiag) / (dtau)
  rescale_cnt = rescale_cnt + 1

lyap_avg = lyap_sum / rescale_cnt
sv.setLEs(lyap_avg)

print(sv)
print('Lyapunov Exponents:')
print(lyap_avg)

#------------------------------------------------------------------
# Store the nature run data
#------------------------------------------------------------------
outfile = name+'_LEs.pkl'
sv.save(outfile)
