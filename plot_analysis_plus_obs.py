from class_lorenz63 import lorenz63
from class_state_vector import state_vector
from class_obs_data import obs_data

#------------------------------------------------------------------
# Read state vector object
#------------------------------------------------------------------
#method = 'skip'
#method = 'nudging'
#method='OI'
method = '3DVar'
infile='x_analysis_'+method+'.pkl'
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
l63 = lorenz63()
error = obs.getVal()-sv.getTrajectory()
#error = abs(obs.getVal()-sv.getTrajectory())
states = sv.getTrajectory()
print('states = ')
print(sv.getTrajectory())
l63.plot_lines_and_points(states=states,points=obs.getVal(),cvec=error)
