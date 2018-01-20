import numpy as np
from class_lorenz63 import lorenz63
from class_state_vector import state_vector
from class_obs_data import obs_data
from class_da_system import da_system

#-----------------------------------------------------------------------
# Read the L63 nature run
#-----------------------------------------------------------------------
infile = 'x_nature.pkl'
sv = state_vector()
sv = sv.load(infile)
x_nature = sv.getTrajectory()

#-----------------------------------------------------------------------
# Read the L63 observations
#-----------------------------------------------------------------------
infile = 'y_obs.pkl'
obs = obs_data()
obs = obs.load(infile)

#-----------------------------------------------------------------------
# Try reducing the observed dimensions
#-----------------------------------------------------------------------
#obs.reduceDim([0])  # x only
#obs.reduceDim([1])  # y only
#obs.reduceDim([2])  # z only
#obs.reduceDim([0,1])  # x and y only
#obs.reduceDim([1,2])  # y and z only
#obs.reduceDim([0,2])  # z and x only

y_obs = obs.getVal()
y_pts = obs.getPos()
y_err = obs.getErr()
print('y_obs = ')
print(y_obs[0,:])
print('y_pts = ')
print(y_pts[0,:])

#-----------------------------------------------------------------------
# Initialize the da system
#-----------------------------------------------------------------------
das = da_system(state_vector = sv, obs_data = obs)

I = np.identity(xdim)

# Set background error covariance
sigma_b = 0.9

# Set observation error covariance
sigma_r = 1.0

# Set the linear observation operator matrix as the identity by default 
H = I

print('B = ')
print(B)
print('R = ')
print(R)
print('H = ')
print(H)

das.setB(sigma_b**2*I)
das.setR(sigma_r**2*I)
das.setH(I)

#-----------------------------------------------------------------------
# Initialize the timesteps
#-----------------------------------------------------------------------
t_nature = sv.getTimes()
ainc_step = 1  # (how frequently to perform an analysis)
dtau = (t_nature[ainc_step] - t_nature[0])
tsteps=10 * ainc_step
dt = dtau/tsteps
maxit,xdim = np.shape(x_nature)

# Store basic timing info in das object
das.ainc = ainc_step
das.dtau = dtau
das.tstep = tsteps
das.dt = dt
das.maxit = maxit
das.xdim = xdim

#-----------------------------------------------------------------------
# Choose DA method:
#-----------------------------------------------------------------------

#-----------
# Test basic functionality
#-----------
method='skip'

#-----------
# 3D methods
#-----------
# Nudging
#method='nudging'
# OI
#method='OI'
# 3D-Var
#method='3DVar'

das.setMethod(method)

#-----------------------------------------------------------------------
# Store DA object
#-----------------------------------------------------------------------
name = 'x_analysis_init'
outfile=name+'.pkl'
das.save(outfile)
