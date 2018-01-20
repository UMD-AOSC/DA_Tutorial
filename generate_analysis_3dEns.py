# Tutorial: "A tour of Data Assimilation methods"
# Model: Lorenz-63
# DA Methods: Nudging, 3D-Var, 4D-Var, Particle Filter, EnKF, Hybrid
import numpy as np
from class_lorenz63 import lorenz63
from class_state_vector import state_vector
from class_obs_data import obs_data
from class_da_system import da_system

#-----------------------------------------------------------------------
# Read the da system object
#-----------------------------------------------------------------------
name = 'x_analysis_init'
infile = name+'.pkl'
sv = das.getStateVector()
x_nature = sv.getTrajectory()

#-----------------------------------------------------------------------
# Initialize the timesteps
#-----------------------------------------------------------------------
t_nature = sv.getTimes()
ainc_step = das.ainc  # (how frequently to perform an analysis)
dtau = das.dtau
tsteps= das.tsteps
dt = das.dt
maxit = das.maxit
xdim = das.xdim

#-----------------------------------------------------------------------
# Get the L63 observations via the obs_data object
#-----------------------------------------------------------------------
obs = das.getObsData()

#-----------------------------------------------------------------------
# Initialize the model
#-----------------------------------------------------------------------
l63 = lorenz63()

#-----------------------------------------------------------------------
# Choose DA method:
#-----------------------------------------------------------------------

method = das.getMethod()  # (use default)

#-----------
# Test basic functionality
#-----------
#method='skip'

#-----------
# Ensemble methods
#-----------
# Particle filter
#method='PF'
# EnKF
#method='ETKF'

das.setMethod(method)

#-----------------------------------------------------------------------
# Initialize the ensemble
#-----------------------------------------------------------------------
xa = sv.x0
bias_init = 0
sigma_init = 0.1
edim = 3 #20
Xa = das.initEns(xa,mu=bias_init,sigma=sigma_init,edim=edim)

#-----------------------------------------------------------------------
# Conduct data assimilation process
#-----------------------------------------------------------------------
#
xa = sv.x0
xa_history = np.zeros_like(x_nature)
KH_history = []
for i in range(0,maxit,ainc_step):
 
  #----------------------------------------------
  # Run forecast ensemble for this analysis cycle:
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
  yp = y_pts[i,:]

  #----------------------------------------------
  # Update the error covariances
  #----------------------------------------------
# das.setB(sigma_b**2*I)
# das.setR(sigma_r**2*I)
# das.setH(I)
# das.reduceYdim(yp)
 
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
  xa_history[i,:] = xa

  # Archive the KH matrix
  KH_history.append(KH)
 
#--------------------------------------------------------------------------------
# Fill in unobserved dimensions (for plotting)
#--------------------------------------------------------------------------------
#fillValue=0.0
#obs.fillDim(fillValue)
#das.setObsData(obs)

sv.setTrajectory(xa_history)
das.setStateVector(sv)

outfile='x_analysis_'+method+'.pkl'
das.save(outfile)
