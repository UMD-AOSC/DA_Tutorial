# Tutorial: "A tour of Data Assimilation methods"
# Model: Lorenz-63
# DA Methods: Nudging, 3D-Var, 4D-Var, Particle Filter, EnKF, Hybrid
import numpy as np
from class_lorenz63 import lorenz63
from class_state_vector import state_vector
from class_obs_data import obs_data
from class_da_system import da_system

#-----------------------------------------------------------------------
# Read the L63 nature run
#-----------------------------------------------------------------------
name = 'x_nature'
infile = name+'_QR.pkl'
sv = state_vector()
sv = sv.load(infile)
x_nature = sv.getTrajectory()
t_nature = sv.getTimes()

#-----------------------------------------------------------------------
# Initialize the timesteps
#-----------------------------------------------------------------------
t_nature = sv.getTimes()
dt = (t_nature[1] - t_nature[0])
ainc_step = sv.hist_ainc
dtau = sv.hist_dtau
_,xdim = np.shape(x_nature)

#------------------------------------------------------------------
# Compute the LEs by propagating the TLM in time and intermittently 
# performing a QR decomposition and rescaling as necessary.
#------------------------------------------------------------------
I = np.identity(xdim)
Mhist = sv.getMhist()
Rhist = sv.getRhist()

lyap_sum = np.zeros(xdim)
rescale_cnt = 0
maxit = len(Rhist)
for i in range(maxit):
  
  # In order to find the eigenvalues of the final matrix describing
  # the evolution of the entire trajectory, we would like to 
  # compute the product of the M matrices to identify the
  # stretching and contracting that results.
  #

  # Interlude: (optional)
  # Try computing the product of all M matrices to see what happens 

  # In order to compute the eigenvalues, we will need to transform
  # and work with summing exponents instead of multiplying matrices
  # The final summation of the exponents, scaled by the time duration
  # of the trajectory, gives the 'Lyapunov Exponents'.
  # The eigenvalues of the matrices M = QR are contained in the 
  # diagonal of the upper-triangular R matrices.
  #
  R = Rhist[i]
  rdiag = abs(np.diag(R))
  lyap_sum = lyap_sum + np.log(rdiag) / (dtau)
  rescale_cnt = rescale_cnt + 1

  # Check in comparison to published results:
  # http://sprott.physics.wisc.edu/chaos/lorenzle.htm

lyap_avg = lyap_sum / rescale_cnt
sv.setLEs(lyap_avg)
print(sv)
print('Lyapunov Exponents:')
print(lyap_avg)

print('The sum of the LE\'s should equal the trace of the Jacobian (−13.6666).')
print('Sum(LEs) = ', np.sum(lyap_avg))

#------------------------------------------------------------------
# Store the nature run data
#------------------------------------------------------------------
outfile = name+'_LEs.pkl'
sv.save(outfile)
