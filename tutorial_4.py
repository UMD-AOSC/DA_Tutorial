# Tutorial: "A tour of Data Assimilation methods"
# Model: Lorenz-63
# DA Methods: Nudging, 3D-Var, 4D-Var, Particle Filter, EnKF, Hybrid
import numpy as np
from class_lorenz63 import lorenz63
from class_state_vector import state_vector
from class_obs_data import obs_data
from class_da_system import da_system

#-----------------------------------------------------------------------
# Exercises:
# (1) Adjust parameters to 'break' the methods:
#  (a) Increase observational noise (by increasing sigma_r)
#  (b) Observe fewer dimensions (e.g. only x and y, only y, only z) by modifying the H operator
#  (b) Add a bias to the model
#  (c) Use a fundamentally different model as 'truth' (i.e. introduce a systematic model error)
#  (d) Draw observational errors from a skewed distribution
# (2) Test the Hybrid filter
#  (a) Run with many and few members
#  (b) Change the definition of the observation error
#  (c) Change the forecast length
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
# Read the L63 nature run
#-----------------------------------------------------------------------
infile = 'x_nature.pkl'
sv = state_vector()
sv = sv.load(infile)
x_nature = sv.getTrajectory()

#-----------------------------------------------------------------------
# Initialize the timesteps
#-----------------------------------------------------------------------
t_nature = sv.getTimes()
ainc_step = 1  # (how frequently to perform an analysis)
dtau = (t_nature[ainc_step] - t_nature[0])
tsteps=10 * ainc_step
dt = dtau/tsteps
maxit,xdim = np.shape(x_nature)

#-----------------------------------------------------------------------
# Read the L63 observations
#-----------------------------------------------------------------------
infile = 'y_obs.pkl'
obs = obs_data()
obs = obs.load(infile)
y_obs = obs.getVal()
y_pts = obs.getPos()
y_err = obs.getErr()
print('y_obs = ')
print(y_obs[0,:])
print('y_pts = ')
print(y_pts[0,:])

#-----------------------------------------------------------------------
# Initialize the model
#-----------------------------------------------------------------------
l63 = lorenz63()

#-----------------------------------------------------------------------
# Initialize the da system
#-----------------------------------------------------------------------
das = da_system()
I = np.identity(xdim)

sigma_b = 0.9
das.setB(sigma_b**2*I)
print('B = ')
print(das.B)

sigma_r = 1.0
das.setR(sigma_r**2*I)
print('R = ')
print(das.R)

das.setH(I)
print('H = ')
print(das.H)

#-----------------------------------------------------------------------
# Choose DA method:
#-----------------------------------------------------------------------

#-----------
# Session 0:
# Test basic functionality
#-----------
#method='skip'

#-----------
# Session 4:
# Hybrid methods
#-----------
# Hybrid gain
method='Hybrid'

das.setMethod(method)

#-----------------------------------------------------------------------
# Initialize the ensemble
#-----------------------------------------------------------------------
xa = sv.x0
bias_init = 0
sigma_init = 0.1
edim = 20
Xa = das.initEns(xa,mu=bias_init,sigma=sigma_init,edim=edim)

#-----------------------------------------------------------------------
# Test assimilation methods:
#-----------------------------------------------------------------------
xa_history = np.zeros_like(x_nature)
KH_history = []
for i in range(0,maxit,ainc_step):
 
  #----------------------------------------------
  # Run forecast model for this analysis cycle:
  #----------------------------------------------
  t = np.arange(t_nature[i],t_nature[i]+dtau+dt,dt)
# print('t = ', t)
  # Run the model ensemble forecast
  Xf = np.zeros_like(Xa)
  for k in range(edim):
    xf_4D =  l63.run(Xa[:,k].A1,t) 
    # Get last timestep of the forecast
    Xf[:,k] = np.transpose(np.matrix(xf_4D[-1,:]))

  #----------------------------------------------
  # Get the observations for this analysis cycle
  #----------------------------------------------
  yo = y_obs[i,:]

  #----------------------------------------------
  # Update the error covariances
  #----------------------------------------------
# das.setB(sigma_b**2*I)
# das.setR(sigma_r**2*I)
# das.setH(I)
 
  #----------------------------------------------
  # Compute analysis
  #----------------------------------------------
  Xa, KH = das.compute_analysis(Xf,yo)
  xa = np.mean(Xa,axis=1)

# print('Xa = ')
# print(Xa)
# print('xa = ')
# print(xa)
# print('xn = ')
# print(x_nature[i,:].T)
# print('KH = ')
# print(KH)

  # Archive the analysis
  xa_history[i,:] = xa.T

  # Archive the KH matrix
  KH_history.append(KH)

# if ((i+1)%5==0):
#   exit()
 
sv.setTrajectory(xa_history)
outfile='x_analysis_'+method+'.pkl'
sv.save(outfile)

