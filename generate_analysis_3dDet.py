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
# Session 0:
# Test basic functionality
#-----------
#method='skip'

#-----------
# Session 2:
# 3D methods
#-----------
# Nudging
#method='nudging'
# OI
method='OI'
# 3D-Var
#method='3DVar'

das.setMethod(method)

#-----------------------------------------------------------------------
# Conduct data assimilation process
#-----------------------------------------------------------------------
#
xa = sv.x0
xa_history = np.zeros_like(x_nature)
KH_history = []
for i in range(0,maxit,ainc_step):
 
  #----------------------------------------------
  # Run forecast model for this analysis cycle:
  #----------------------------------------------
  t = np.arange(t_nature[i],t_nature[i]+dtau+dt,dt)
# print('t = ', t)
  # Run the model
  xf_4D =  l63.run(xa,t) 
  # Get last timestep of the forecast
  xf = xf_4D[-1,:] 

  #----------------------------------------------
  # Get the observations for this analysis cycle
  #----------------------------------------------
  yo = y_obs[i,:]
  yp = y_pts[i,:]

  #----------------------------------------------
  # Update the error covariances
  #----------------------------------------------
  das.setB(sigma_b**2*I)
  das.setR(sigma_r**2*I)
  das.setH(I)
  das.reduceYdim(yp)
 
  #----------------------------------------------
  # Compute analysis
  #----------------------------------------------
  xa, KH = das.compute_analysis(xf,yo)
# print('xa = ')
# print(xa)
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
