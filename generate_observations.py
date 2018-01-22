import numpy as np
from class_state_vector import state_vector
from class_obs_data import obs_data

infile = 'x_nature.pkl'
outfile = 'y_obs.pkl'

mu = 0
sigma = 0.001

#--------------------------------------------------------------------------------
# Create observation object
#--------------------------------------------------------------------------------
obs = obs_data(name='observe_full_state',err=[],bias=mu,val=[],pos=[])

#--------------------------------------------------------------------------------
# Read the nature run
#--------------------------------------------------------------------------------
sv = state_vector()
sv = sv.load(infile)
x_nature = sv.getTrajectory()

nr,nc = x_nature.shape
print('nr = ',nr)
print('nc = ',nc)

#--------------------------------------------------------------------------------
# Compute the climatological variance
#--------------------------------------------------------------------------------
x_std = np.zeros(nc)
for i in range(nc):
  x_std[i] = np.std(x_nature[:,i])

#--------------------------------------------------------------------------------
# Sample the nature run and apply noise
#--------------------------------------------------------------------------------
yo = np.zeros_like(x_nature)
eta = np.zeros_like(x_nature)
hx = np.zeros_like(x_nature)
for i in range(nc):
  # Compute error as a percentage of climatological variance
  eta[:,i] = np.random.normal(mu,sigma*x_std[i],nr)
  # (Could apply H(x_nature) here):
  hx[:,i] = x_nature[:,i]  
  # Observation is simulated as H(x_nature) + noise
  yo[:,i] = hx[:,i] + eta[:,i]

pos = np.zeros_like(yo)
for i in range(nr):
  pos[i,:] = [0,1,2]

obs.setVal(yo)
obs.setErr(eta)
obs.setHx(hx)
obs.setDep(yo-hx)
obs.setPos(pos)

#--------------------------------------------------------------------------------
# Store the true and noisy observations
#--------------------------------------------------------------------------------
obs.save(outfile)

