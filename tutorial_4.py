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
# (1) Adjust parameters to 'break' the methods:
#  (a) Increase observational noise (by increasing sigma_r)
#  (b) Observe fewer dimensions (e.g. only x and y, only y, only z) by modifying the H operator
#  (b) Add a bias to the model
#  (c) Use a fundamentally different model as 'truth' (i.e. introduce a systematic model error)
#  (d) Draw observational errors from a skewed distribution
# (2) Test the Hybrid filter
#  (a) Run with many and few members
#  (b) Change the definition of the observation error
#  (c) Change the forecast length
#-----------------------------------------------------------------------

