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
# (1) Compute the Lyapunov exponents of the nature system
# (2) Compute the forward and backward Lyapunov vectors
# (3) Compute the Covariant Lyapunov Vectors
#-----------------------------------------------------------------------


#-----------------------------------------------------------------------
# Read the L63 nature run
#-----------------------------------------------------------------------
name = 'x_nature'
infile = name+'.pkl'
sv = state_vector()
sv = sv.load(infile)
x_nature = sv.getTrajectory()
t_nature = sv.getTimes()
tsteps=10
dtau = (t_nature[1] - t_nature[0])
dt = dtau/tsteps
maxit,xdim = np.shape(x_nature)
sigma,rho,beta = sv.params

#------------------------------------------------------------------
# Initialize the l63 object
#------------------------------------------------------------------
l63 = lorenz63(sigma=sigma,rho=rho,beta=beta)

#------------------------------------------------------------------
# Input a nature run and compute its corresponding Jacobian
#------------------------------------------------------------------
print('Compute Jacobian at each timestep...')
Jhist = l63.compute_Jfd(sv.getTrajectory(),sv.getTimes())
sv.setJhist(Jhist)

#------------------------------------------------------------------
# Store the nature run with corresponding history of the Jacobian
#------------------------------------------------------------------
outfile = name+'_Jhist.pkl'
sv.save(outfile)
