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
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
# Tutorial 2:
#
# Exercises:
# (1) Generate an analysis using nudging
# (2) Generate an analysis using Optimal Interpolation (OI)
# (3) Generate an analysis using 3D-Var
#
# Additional exercises:
# (1) Test different estimates of the background error covariance matrix
# (2) Test different estimates of the observation error
# (3) Adjust parameters to 'break' the methods:
#  (a) Increase observational noise (by increasing sigma_r)
#  (b) Observe fewer dimensions (e.g. only x and y, only y, only z) by modifying the H operator
#  (c) Add a bias to the model
#  (d) Use a fundamentally different model as 'truth' (i.e. introduce a systematic model error)
#  (e) Draw observational errors from a skewed distribution
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
# Step 2a:
#-----------------------------------------------------------------------
# Run the python program:
# python generate_analysis_3dDet.py
#
# This program generates an analysis using a 3D deterministic method.
# The options for DA method include:
#   'skip'    :: skips the DA procedure
#   'nudging; :: applies a procedure directly forcing the state 
#                to the observations
#   'OI'      :: Optimal interpolation
#   '3D-Var'  :: 3D-Variational method
#   
# Start with option method='skip', plot the results. It should match
# the nature run unless you add a small perturbation to the initial
# conditions. Try different sized perturbations to the initial conditions.
#


#-----------------------------------------------------------------------
# Step 2b:
#-----------------------------------------------------------------------
# Run the python program:
# python generate_analysis_3dDet.py
#
# Try method='nudging', plot the results with:
#
# Try varying:
#  Initial perturbation
#  Observation error (sigma_r)
#  Magnitude of constant nudging factor 
#


#-----------------------------------------------------------------------
# Step 2c:
#-----------------------------------------------------------------------
# Run the python program:
# python generate_analysis_3dDet.py
#
# Try method='OI', plot the results with:
#
# Try varying:
#  Initial perturbation
#  Observation error (sigma_r)
#  Magnitude of constant nudging factor 
#


#-----------------------------------------------------------------------
# Step 2d:
#-----------------------------------------------------------------------
# Run the python program:
# python generate_analysis_3dDet.py
#
# Try method='3DVar', plot the results with:
#
# Try varying:
#  Initial perturbation
#  Observation error (sigma_r)
#  Magnitude of constant nudging factor 
#


#-----------------------------------------------------------------------
# Step 3:
#-----------------------------------------------------------------------
#
# Compare the methods 'OI' and '3DVar' with identical input parameters.
#
# Are the results the same? Should they be?
#
