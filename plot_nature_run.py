from class_lorenz63 import lorenz63
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
l63 = lorenz63()
l63.plot(sv.getTrajectory(),sv.t)
