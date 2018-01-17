import numpy as np
from class_lorenz63 import lorenz63
from class_state_vector import state_vector

#-----------------------------------------------------------------------
# Read the L63 nature run to use identical parameters
#-----------------------------------------------------------------------
infile = 'x_nature.pkl'
sv = state_vector()
sv = sv.load(infile)
x_nature = sv.getTrajectory()
t_nature = sv.getTimes()
maxit,xdim = np.shape(x_nature)
sigma,rho,beta = sv.params

#------------------------------------------------------------------
# Initialize the l63 object
#------------------------------------------------------------------
l63 = lorenz63(sigma=sigma,rho=rho,beta=beta)

#------------------------------------------------------------------
# Run L63 to generate a nature run and its corresponding TLM
#------------------------------------------------------------------
print('Run l63 and compute TLM at each timestep...')
print('sv = ')
print(sv)
Jhist = l63.compute_Jfd(sv.trajectory,sv.t)

# Store the history of TLMs
sv.setTLM(Jhist)

#------------------------------------------------------------------
# Store the nature run data
#------------------------------------------------------------------
outfile = 'x_nature_TLM.pkl'
sv.save(outfile)
