import numpy as np
from scipy.integrate import odeint
#import plotly.offline as py
#import plotly.graph_objs as go
#from plotly import tools
from copy import deepcopy

#===============================================================================
# Define functions used by the class
#===============================================================================

#-------------------------------------------------------------------------------
def f(x, t, F=8.0):
#-------------------------------------------------------------------------------
  # x is the state vector, of any length >= 3
  # F is the forcing, which can either be a scalar or a vector of size equal to len(x)

  dN = len(x)
  # Compute state derivatives
  dx = np.zeros(dN)

  # Compute the 3 edge cases: i=1,2,N
  dx[   0] = (x[1] - x[dN-2]) * x[dN-1] - x[  0]
  dx[   1] = (x[2] - x[dN-1]) * x[   0] - x[  1]
  dx[dN-1] = (x[0] - x[dN-3]) * x[dN-2] - x[dN-1]
  
  # Then the general case
#  for i in range(2, D-1):
#    dx[i] = (x[i+1] - x[i-2]) * x[i-1] - x[  i]

  # Use slice function for the rest of the vector to facilitate improved scaling:
  si = slice(2,dN-1)
  sip1 = slice(3,dN)
  sim2 = slice(0,dN-3)
  sim1 = slice(1,dN-2)
    
  dx[si] = (x[sip1]-x[sim2])*x[sim1]-x[si]

  # Add the forcing
  dx = dx + F

  return dx  
  

#-------------------------------------------------------------------------------
def Ja(x, t, params):
#-------------------------------------------------------------------------------
  # Compute the analytic Jacobian of the Lorenz-96 system
  # at a point represented by 'state'. The time 't' is unused.
  F=params[0]
  dN = len(x)
  # Set up -1.0 on the diagonal
  J = -np.eye(dN)

  # Construct each row
  for i in range(dN):
    jm2 = np.mod(i-2,dN)
    jm1 = np.mod(i-1,dN)
    jp1 = np.mod(i+1,dN)
    J[i,jm2] = -x[jm1]
    J[i,jm1] = x[jp1] - x[jm1]
    J[i,jp1] = x[jm1] 

  J = np.matrix(J)
  return J

#===============================================================================
class lorenz96:
#===============================================================================

  #-----------------------------------------------------------------------------
  # Initialize model parameters
  #-----------------------------------------------------------------------------
  def __init__(self, F=8.0):
    self.F = F
    self.params = [self.F]

  #-----------------------------------------------------------------------------
  # Run model
  #-----------------------------------------------------------------------------
  def run(self, state0, t):
    states = odeint(f, state0, t, args=(self.F,))
    return states

  #-----------------------------------------------------------------------------
  # Compute approximate TLM with I + Df(x0)*dt
  #-----------------------------------------------------------------------------
  def compute_TLMa(self, states, t):

#   print('states = ')
#   print(states)

#   print('times = ')
#   print(t)

    nr,nc = np.shape(states)
    I = np.identity(nc)

    # Compute linear propagator for each timestep
    maxit = len(t)
    Mhist=[]
    for i in range(maxit):
      if (i < maxit-1):
        dt = t[i+1] - t[i]
      else:
        dt = t[-1] - t[-2]

      # Evaluate Jacobian
      Df = Ja(states[i,:], t[i], self.params)

#     print('Df = ')
#     print(Df)

      # Compute approximate linear propagator
      # (Note this is a poor approximation of the matrix integral, 
      # we would prefer a geometric integrator, e.g. the Magnus expansion)
      M = I + Df*dt  
      Mhist.append(deepcopy(M))

    return Mhist

  #------------------------------------------------------------------
  # Plot model output
  #------------------------------------------------------------------
  def plot(self, states,cvec,outfile='l96-3d',plot_title='Lorenz 96 attractor'):
    pass

  #-----------------------------------------------------------------------------
  # Plot model output lines and points
  #-----------------------------------------------------------------------------
  def plot_lines_and_points(self,states,points,cvec,outfile='l96-3d',plot_title='Lorenz 96 attractor'):
    pass

  #-----------------------------------------------------------------------------
  # Plot model output lines and lines
  #-----------------------------------------------------------------------------
  def plot_lines_and_lines(self,states1,states2,cvec,outfile='l96-3d',plot_title='Lorenz 96 attractor',name1='trajectory 1', name2='trajectory 2'):
    pass


