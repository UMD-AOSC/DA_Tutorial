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

#-----------------------------------------------------------------------
# Step 1:
#-----------------------------------------------------------------------
# Run the python program:
# python analysis_init.py Hybrid
# python generate_analysis_3dEns.py
#
# Try method='Hybrid', plot the results.
#
# Vary the parameters of the hybrid and compare to previous methods.
# Try changing the forecast length (acyc_step in analysis_init.py),
# the ensemble size, the number of observations assimilated, the amount
# of model bias or parameter error.
#
# python plot_analysis_vs_analysis.py Hybrid ETKF
# python plot_analysis_vs_analysis.py Hybrid 3DVar
# python plot_analysis_vs_analysis.py Hybrid PF
#
