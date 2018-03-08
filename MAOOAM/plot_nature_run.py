from class_maooam import maooam
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
model = maooam()
xidx = [0,7,14]
model.plot(sv.getTrajectory(),sv.t,xidx=xidx)
