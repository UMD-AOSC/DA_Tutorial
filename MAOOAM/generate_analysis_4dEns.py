# Tutorial: "A tour of Data Assimilation methods"
# Model: MAOOAM
# DA Methods: Nudging, 3D-Var, 4D-Var, Particle Filter, EnKF, Hybrid
import numpy as np
from class_maooam import maooam
from class_state_vector import state_vector
from class_obs_data import obs_data
from class_da_system import da_system
from copy import deepcopy

#-----------------------------------------------------------------------
# Read the da system object
#-----------------------------------------------------------------------
name = 'x_analysis'
infile = name+'_init.pkl'
das = da_system()
das = das.load(infile)

print(das)

#-----------------------------------------------------------------------
# Initialize the ensemble
#-----------------------------------------------------------------------
edim = das.edim
Xa = das.initEns(das.x0,mu=das.ens_bias_init,sigma=das.ens_sigma_init,edim=das.edim)

print('ensemble dimension = ')
print(das.edim)
print('initial bias = ')
print(das.ens_bias_init)
print('initial standard deviation = ')
print(das.ens_sigma_init)
print('X0 = ')
print(Xa)

#-----------------------------------------------------------------------
# Get the nature run trajectory
#-----------------------------------------------------------------------
sv = das.getStateVector()
x_nature = sv.getTrajectory()

#-----------------------------------------------------------------------
# Get the MAOOAM observations via the obs_data object
#-----------------------------------------------------------------------
obs = das.getObsData()
y_obs = obs.getVal()
y_pts = obs.getPos()
y_err = obs.getErr()
print('y_obs = ')
print(y_obs[0,:])
print('y_pts = ')
print(y_pts[0,:])

#-----------------------------------------------------------------------
# Initialize the timesteps
#-----------------------------------------------------------------------
t_nature = sv.getTimes()
acyc_step = das.acyc_step  # (how frequently to perform an analysis)
dtau = das.dtau
dt = das.dt
fcst_step= das.fcst_step
fcst_dt = das.fcst_dt
maxit = das.maxit
xdim = das.xdim
ydim = das.ydim

#-----------------------------------------------------------------------
# Initialize the model
#-----------------------------------------------------------------------
model = maooam()

#-----------------------------------------------------------------------
# Choose DA method:
#-----------------------------------------------------------------------
method = das.getMethod()

#-----------------------------------------------------------------------
# Conduct data assimilation process
#-----------------------------------------------------------------------
#
xa = das.x0
xa_history = np.zeros_like(x_nature)
xa_history[:] = np.nan
KH_history = []
KH_idx = []
for i in range(0,maxit-acyc_step,acyc_step):
 
  #----------------------------------------------
  # Get the observations for this analysis cycle
  #----------------------------------------------
  yo_4d = y_obs[i:i+acyc_step+1,:]
  yp_4d = y_pts[i:i+acyc_step+1,:]
  tdim,ydim = np.shape(yo_4d)

  #----------------------------------------------
  # Run forecast model for this analysis cycle:
  #----------------------------------------------
  t = np.linspace(t_nature[i],t_nature[i+acyc_step], acyc_step+1, endpoint=True)
# print('t = ', t)
  if (len(t)>acyc_step+1):
    print('acyc_step = ', acyc_step)
    print('dt = ', dt)
    print('t_nature[i] = ', t_nature[i])
    print('t_nature[i+acyc_step] = ', t_nature[i+acyc_step])
    print('t_nature[i+acyc_step+1] = ' ,t_nature[i+acyc_step+1])
    print('len(t) = ',len(t))
    exit()

  #----------------------------------------------
  # Run the model ensemble forecast
  #----------------------------------------------
  X  = np.zeros_like(Xa)
  Xf_4d = []
  # Initialize the 4D list of Xf ensemble forecast matrices
  for ti in range(tdim):
    Xf_4d.append(X)
    
  xf_4d = 0
  # Preferably, run this loop in parallel:
  for k in range(das.edim):
    # Run model run for ensemble member k
    xf_4d_k =  model.run(Xa[:,k].A1,t) 
    print('k = ', k)
    print('xf_4d_k = ')
    print(xf_4d_k)
    # Get each timestep of the forecast
    for ti in range(tdim):
      print('ti = ', ti)
      print('Xf_4d = ')
      print(Xf_4d)
      Xf_4d[ti][:,k] = np.transpose(np.matrix(xf_4d_k[ti,:]))
      print('Xf_4d[ti] = ')
      print(Xf_4d[ti])
      print('Xf_4d[ti][:,k] = ')
      print(Xf_4d[ti][:,k])
    # Compute forecast ensemble mean
    xf_4d = xf_4d + xf_4d_k
  xf_4d = xf_4d / das.edim 

  exit()
  
# print('xf_4d = ')
# print(xf_4d)

  #----------------------------------------------
  # Construct the M,H,R matrix lists.
  # Could potentially have a time-dependent
  # H and R matrix here if desired.
  #----------------------------------------------
# Mlist = model.compute_TLMa(xf_4d,t)
# print('Mlist = ')
# print(Mlist)
  Hlist = das.expandToList(das.H,len(t))
  Rlist = das.expandToList(das.R,len(t))
  params = [Hlist,Rlist]

  #----------------------------------------------
  # Compute analysis
  #----------------------------------------------
  Xa, KH = das.compute_analysis(Xf_4d,yo_4d,params)
  xa = np.mean(Xa,axis=1).T

# print('xa = ')
# print(xa)
# print('x_nature = ')
# print(x_nature[i+acyc_step,:])
# print('KH = ')
# print(KH)

# exit('exiting on purpose... post-analysis.')

  # Fill in the missing timesteps with the forecast from the previous analysis IC's
  xa_history[i:i+acyc_step,:] = xf_4d[0:acyc_step,:]
  # Archive the analysis
  xa_history[i+acyc_step,:] = xa

# print('xa_history[i:i+acyc_step+1,:] = ', xa_history[i:i+acyc_step+1,:])

  # Archive the KH matrix
  KH_history.append(deepcopy(KH))
  KH_idx.append(i+acyc_step)
 
das.setKH(KH_history,KH_idx)

print('Last background error covariance matrix B = ')
print(das.getB())

sv.setTrajectory(xa_history)
sv.setName(name)
das.setStateVector(sv)

outfile=name+'_'+method+'.pkl'
das.save(outfile)

print(das)
