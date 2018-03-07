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
# (1) Test the EnKF
#  (a) Run with many and few members
#  (b) Change the definition of the observation error
#  (c) Change the forecast length
# (2) Test the particle filter
#  (a) Run with many and few members
#  (b) Change the definition of the observation error
#  (c) Change the forecast length
#
# Additional Exercises:
# (1) Adjust parameters to 'break' the methods:
#  (a) Increase observational noise (by increasing sigma_r)
#  (b) Observe fewer dimensions (e.g. only x and y, only y, only z) by modifying the H operator
#  (c) Add a bias to the model
#  (d) Use a fundamentally different model as 'truth' (i.e. introduce a systematic model error)
#  (e) Draw observational errors from a skewed distribution
# (2) Explore the use of 'inflation'
#-----------------------------------------------------------------------


#-----------------------------------------------------------------------
# Step 1:
#-----------------------------------------------------------------------
# Run the python program:
# python generate_analysis_init.py
#
# Set up the initial conditions as desired for the data assimilation
# experiment. Background and observation error covariance matrices,
# and observation operator can be modified here. Also, analysis increment
# duration, and model foreast dt.
#


#-----------------------------------------------------------------------
# Step 2:
#-----------------------------------------------------------------------
# Run the python program:
# python generate_analysis_3dEns.py
#
# This program generates an analysis using a 3D Ensemble method.
# The options for DA method include:
#   'skip'    :: skips the DA procedure
#   'PF'      :: applies a simple Sequential/Sampling Importance 
#                Resampling (SIR) particle filter
#   'ETKF'    :: applies Ensemble Transform Kalman Filter
#   
# Start with option method='skip', plot the results. Does it match the 
# nature run? Should it?
#
# Experiment with different initial perturbations to generate the initial
# ensemble. Try adding a bias to the initial ensemble.
#


#-----------------------------------------------------------------------
# Step 3:
#-----------------------------------------------------------------------
# Run the python program:
# python generate_analysis_3dEns.py
#
# Use method='ETKF', plot the results.
#
# As before, try different initial conditions for the ensemble.
#
# The ETKF is ideal when the distribution can be approximated as
# Gaussian, usually when the timescale of the analysis cycle exhibits
# roughly linear error dynamics.
#
# The ETKF automatically estimates the background error covariance as
# a time varying estimate B = (1/n)*(Xb * Xb.T), where Xb is the set of 
# perturbations from the ensemble mean and 'n' is the ensemble size. 
# Under what conditions might this approximation fail? Why?
#
# Try using various ensemble sizes and analysis cycle lengths. Also, try
# adding a bias to the forecast model. Can you find conditions that 
# break the ETKF?
#
# Try increasing the observation error. Recall that the 'true' error
# does not have to be identical to the estimate error (R). Is it better
# to overestimate or underestimate the observation error?
#


#-----------------------------------------------------------------------
# Step 4:
#-----------------------------------------------------------------------
# Run the python program:
# python generate_analysis_3dEns.py
#
# Use method='PF', plot the results.
#
# As before, try different initial conditions for the ensemble.
#
# The PF is generally considered ideal for highly nonlinear dynamics or
# non-Guassian distributions. Try increasing the analysis cycle length
# to push the limits of the PF as the DA becomes less constrained by obs.
#
# Try using various ensemble sizes.
#
# Try increasing the observation error. Unlike the nudging/OI/Var methods,
# there is no mechanism in the basic SIR PF to push the state towards
# observations (and thus away from the attractor with noisy obs). This
# implementation does use an additive noise to prevent filter collapse.
#
# (Advanced) Try turning off the inflation in the PF implementation.
# Plot the results. What happens?
#
# (Advanced) Try changing the distribution of the observation error. Can
# you break the ETKF but maintain stability with the PF?
#


