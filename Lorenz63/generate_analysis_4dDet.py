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
# Read the da system object
#-----------------------------------------------------------------------
name = 'x_analysis'
infile = name+'_init.pkl'
das = da_system()
das = das.load(infile)

print(das)

#-----------------------------------------------------------------------
# Get the nature run trajectory
#-----------------------------------------------------------------------
sv = das.getStateVector()
x_nature = sv.getTrajectory()

#-----------------------------------------------------------------------
# Get the L63 observations via the obs_data object
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
l63 = lorenz63()

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
  # Run forecast model for this analysis cycle:
  #----------------------------------------------
# t = np.arange(t_nature[i],t_nature[i+acyc_step+1],dt)  #NOTE: gives inconsistent results.
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

  # Run the model
  xf_4d =  l63.run(xa,t) 
  
# print('xf_4d = ')
# print(xf_4d)

  #----------------------------------------------
  # Get the observations for this analysis cycle
  #----------------------------------------------
  yo_4d = y_obs[i:i+acyc_step+1,:]
  yp_4d = y_pts[i:i+acyc_step+1,:]

  #----------------------------------------------
  # Construct the M,H,R matrix lists.
  # Could potentially have a time-dependent
  # H and R matrix here if desired.
  #----------------------------------------------
  Mlist = l63.compute_TLMa(xf_4d,t)
# print('Mlist = ')
# print(Mlist)
  Hlist = das.expandToList(das.H,len(t))
  Rlist = das.expandToList(das.R,len(t))

  #----------------------------------------------
  # Compute analysis
  #----------------------------------------------
  x0_bg = xf_4d[0,:]
  params = [x0_bg, Mlist, Hlist, Rlist]
  for j in range(das.outer_loops):

#   print('pre-analysis outerloop = ', j,' xf_4d = ') 
#   print(xf_4d)

    # Run one outer loop of 4D-Var:
    xr,xc = xf_4d.shape
    yr,yc = yo_4d.shape
    if (xr != yr):
      print ('xr = ', xr)
      print ('yr = ', yr)
      exit('EXITING...')

    xa, KH = das.compute_analysis(xf_4d,yo_4d,params)

#   print('post-analysis outerloop = ', j,' xa = ') 
#   print(xa)

    # Rerun the model
    xf_4d =  l63.run(xa,t) 

  # Set analysis as end of last outer loop forecast
  xa = xf_4d[-1,:]

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
