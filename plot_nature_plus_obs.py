import l63
from class_state_vector import state_vector
from class_obs_data import obs_data

#------------------------------------------------------------------
# Read state vector object
#------------------------------------------------------------------
infile='x_nature.pkl'
sv = state_vector()
sv = sv.load(infile)

print(sv)

#------------------------------------------------------------------
# Read obs data object
#------------------------------------------------------------------
infile='y_obs.pkl'
obs = obs_data()
obs = obs.load(infile)

print(obs)

#------------------------------------------------------------------
# Plot the result
#------------------------------------------------------------------
error = obs.getVal()-sv.getTrajectory()
l63.plot_lines_and_points(states=sv.getTrajectory(),points=obs.getVal(),cvec=error)
