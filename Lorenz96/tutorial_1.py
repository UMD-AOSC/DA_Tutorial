# Tutorial: "A tour of Data Assimilation methods"
# Model: Lorenz-63
# DA Methods: Nudging, 3D-Var, 4D-Var, Particle Filter, EnKF, Hybrid
#
# Purpose:
#  Tutorial 1 focuses on generating the nature run, free run, and
#  Lyapunov Exponents of the nature system, the Lorenz-63 model.
#
#

import numpy as np
from class_lorenz63 import lorenz63
from class_state_vector import state_vector
from class_obs_data import obs_data
from class_da_system import da_system

#-----------------------------------------------------------------------
# Tutorial 1:
#
# Exercises:
# (1) Compute a nature run
# (2) Compute a free run with perturbed initial conditions
# (3) Compute the Jacobian of the system for its entire trajectory
# (4) Compute the Lyapunov exponents of the nature system
# (5) Compute the observations to be used by the data assimilation system
#
# Additional exercises:
# (1) Repeat with different model parameters
# (2) Compute the forward and backward Lyapunov vectors
# (3) Compute the Covariant Lyapunov Vectors
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
# Step 1a:
#-----------------------------------------------------------------------
# Run the python program:
# python generate_nature_run.py
#
# This will generate a file containing the nature run as a 'state_vector'
# object. The object contains the nature run trajectory as well as other
# important information about the nature run and its model parameters,
# similarly to a netcdf file.
#
# Run the python program:
# plot_nature_run.py
#
# Explore the dataset.
#


#-----------------------------------------------------------------------
# Step 1b:
#-----------------------------------------------------------------------
# Generate a free run.
#
# In generate_nature_run.py, Add a small perturbation to the initial 
# conditions to investigate the impacts on the entire model trajectory.
#
# Use the default model parameters, and change the name of the new
# output file from 'x_nature' to 'x_freerun'
# 
# Run the python program:
# python generate_nature_run.py
#
# Run the python program:
# plot_nature_vs_freerun.py
#
# Explore the datasets and compare.
#


#-----------------------------------------------------------------------
# Step 2:
#-----------------------------------------------------------------------
#
# By hand, determine the Jacobian of the Lorenz-63 system.
#
# Run the python program:
# python generate_nature_Mhist.py
#
# This will generate a file containing the same nature run state vector
# with the addition of a history of linear propagator matrices estimated
# by M = I + Df*dt. Different estimates of the linear propagator
# may produce different results due to numerical errors.
#
# Examine the output of the computed M matrices and compare to your
# analytic derivation of the Jacobian matrix.
#


#-----------------------------------------------------------------------
# Step 3:
#-----------------------------------------------------------------------
# Run the python program:
# python generate_nature_QR.py
#
# This will compute a QR decomposition of the linear propagators M = QR
# for a specified time interval over the span of the nature run  trajectory
# (Mhist), and store a history of the Q and R matrices (Qhist, Rhist).
#

#-----------------------------------------------------------------------
# Step 4:
#-----------------------------------------------------------------------
# Run the python program:
# python generate_nature_LEs.py
#
# This will generate the Lyapunov Exponents of the system using the
# information stored in Rhist, the history of R matrices from the QR
# decomposition of the TLM


#-----------------------------------------------------------------------
# Interlude:  (optional)
#-----------------------------------------------------------------------
#
# Further investigation:
# Go back and repeat the process using different model parameters.
# Try a much longer and a much shorter run of the model.
# How do these modifications impact the Lyapunov exponents?
#


#-----------------------------------------------------------------------
# Step 5:
#-----------------------------------------------------------------------
# Run the python program:
# python generate_observations.py
#
# This will generate a set of observations at every point along the 
# trajectory of the nature run by sampling the nature run and applying
# noise scaled relative the the climatological variability of each 
# state element.
#


#-----------------------------------------------------------------------
# Step 6:  (optional)
#-----------------------------------------------------------------------
# Run the python program:
# python generate_nature_CLVs.py
#
# This will generate the Lyapunov Vectors of the system using the
# information stored in Qhist and Rhist, the history of Q and R matrices 
# from the QR decomposition of the TLM
#
