from class_lorenz96 import lorenz96
from class_state_vector import state_vector

#------------------------------------------------------------------
# Read state vector object
#------------------------------------------------------------------
infile='x_nature.pkl'
sv = state_vector()
sv = sv.load(infile)

print(sv)

#------------------------------------------------------------------
# Plot the result
#------------------------------------------------------------------
l96 = lorenz96()
l96.plot(sv.getTrajectory(),sv.t)
