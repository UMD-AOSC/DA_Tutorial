#import numpy as np
from scipy.integrate import odeint
import plotly.offline as py
import plotly.graph_objs as go
from plotly import tools
from copy import deepcopy
from ctypes import *

import numpy as np
import params_maooam
from params_maooam import ndim, tw, t_run, t_trans, dt
#from maooam import params_maooam
#from maooam.params_maooam import ndim, tw, t_run, t_trans, dt
import tl_ad
from maooam import integrator
import time
from maooam import ic_def
from maooam import ic
import sys

class MaooamFortran:
    module_maooam = np.ctypeslib.load_library("step_maooam.so", ".")
    module_maooam.step_maooam_.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64)]
    module_maooam.step_maooam_.restype = c_void_p

    def __init__(self, dt):
        assert dt.__class__ in [float, np.float32, np.float64]
        self.dt = np.array([dt])

    def step(self, x0):
        self.module_maooam.step_maooam_(x0, self.dt)
        return x0

#===============================================================================
# Define functions used by the class
#===============================================================================

# Wrapper for the MAOOAM code to:
# 1) Evolve the evolution function
# 2) Compute the Jacobian
# 3) Compute the corresponding TLM
# 4) Plot slices of the model

#-------------------------------------------------------------------------------
def f(state, t):
#-------------------------------------------------------------------------------
  # Call to MAOOAM tendency function
  # at a point represented by 'state'. The time 't' is unused.
  dxdt = integrator.tendencies(state)
  return dxdt

#-------------------------------------------------------------------------------
def Ja(state, t):
#-------------------------------------------------------------------------------
  # Compute the analytic Jacobian of the MAOOAM system
  # at a point represented by 'state'. The time 't' is unused.
  J = tl_ad.jacobi_mat(state)
  return J

#-------------------------------------------------------------------------------
def Jfd(state0,state1):
#-------------------------------------------------------------------------------
  # Compute the finite difference Jacobian of the MAOOAM system
  # from point 'state0' to point 'state1'

  # If the states are not different, then use machine epsilon as a perturbation
  dist = np.linalg.norm(state1-state0)
  eps = np.finfo(float).eps
  if (dist < eps):
    x_pert = eps*np.ones_like(state1)
    state1 = state0 + x_pert
    print('state0 = ')
    print(state0)
    print('state1 = ')
    print(state1)

  f0 = f(state0,0)
  for j in range(ndim):
    # Initialize at the control state
    xj = state0
    # Perturb only one direction
    xj[j] = state1[j]
    # Evaluate the function at the perturbed state
    fj = f(xj,0)
    # Compute entry of the Jacobian
    dfdx[:,j] = (fj - f0)/(state1[j] - state0[j])

  J  = np.array(dfdx)

  return J


#===============================================================================
class maooam:
#===============================================================================

  #-----------------------------------------------------------------------------
  # Initialize model parameters
  #-----------------------------------------------------------------------------
  def __init__(self, params=[]):
    self.params = params

  #-----------------------------------------------------------------------------
  # Run model
  #-----------------------------------------------------------------------------
  def run(self, state0, t):
#   states = odeint(f, state0, t, args=(self.sigma, self.rho, self.beta))
    # Loop and store a history of states:
    tdim = len(t)
    xdim = len(state0) 
    states = np.zeros((tdim,xdim))
    state = state0
    mf = MaooamFortran(dt)
    for i in range(tdim):
      states[i,:] = state
      state = mf.step(state)

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

    # Compute Jacobian / linear propagator for each timestep
    # sigma,rho,beta = self.params
    maxit = len(t)
    Mhist=[]
    for i in range(maxit):
      if (i < maxit-1):
        dt = t[i+1] - t[i]
      else:
        dt = t[-1] - t[-2]

      # Evaluate Jacobian
      Df = Ja(states[i,:], t[i])

#     print('Df = ')
#     print(Df)

      # Compute approximate linear propagator
      # (Note this is a poor approximation of the matrix integral,
      # we would prefer a geometric integrator, e.g. the Magnus expansion)
      M = I + Df*dt
      Mhist.append(deepcopy(M))

    return Mhist

  #-----------------------------------------------------------------------------
  # Compute approximate TLM with finite differences
  #-----------------------------------------------------------------------------
  def compute_TLMfd(self, states, t):

    print('states = ')
    print(states)

    print('times = ')
    print(t)

    # Compute Jacobian / linear propagator for each timestep
    Jhist=[]
    for i in range(len(t)-1):
      Df = Jfd(states[i,:],states[i+1,:],self.params)

      # Compute approximate linear propagator
      # (Note this is a poor approximation of the matrix integral,
      # we would prefer a geometric integrator, e.g. the Magnus expansion)
      M = I + Df*dt
      Mhist.append(deepcopy(M))

    return Mhist


  #------------------------------------------------------------------
  # Plot 3D-slice of model output
  #------------------------------------------------------------------
  def plot(self, states,cvec,outfile='maooam-3d',plot_title='3D-slice of MAOOAM attractor',xidx=[0,1,2]):

    nr,nc = np.shape(states)
    x = states[:,xidx[0]]
    y = states[:,xidx[1]]
    z = states[:,xidx[2]]

    trace = go.Scatter3d(
      x=x, y=y, z=z,
      marker=dict(
          size=2,
          color=cvec,
          colorscale='Viridis',
      ),
      line=dict(
          color='#1f77b4',
          width=1
      )
    )

    data = [trace]

    layout = dict(
      width=1200,
      height=700,
      autosize=False,
      title=plot_title,
      scene=dict(
          xaxis=dict(
              gridcolor='rgb(255, 255, 255)',
              zerolinecolor='rgb(255, 255, 255)',
              showbackground=True,
              backgroundcolor='rgb(230, 230,230)'
          ),
          yaxis=dict(
              gridcolor='rgb(255, 255, 255)',
              zerolinecolor='rgb(255, 255, 255)',
              showbackground=True,
              backgroundcolor='rgb(230, 230,230)'
          ),
          zaxis=dict(
              gridcolor='rgb(255, 255, 255)',
              zerolinecolor='rgb(255, 255, 255)',
              showbackground=True,
              backgroundcolor='rgb(230, 230,230)'
          ),
          camera=dict(
              up=dict(
                  x=0,
                  y=0,
                  z=1
              ),
              eye=dict(
                  x=-1.7428,
                  y=1.0707,
                  z=0.7100,
              )
          ),
          aspectratio = dict( x=1, y=1, z=1 ),
          aspectmode = 'manual'
      ),
    )

    fig = dict(data=data, layout=layout)
    py.plot(fig, filename=outfile, validate=False)

  #-----------------------------------------------------------------------------
  # Plot model output lines and points
  #-----------------------------------------------------------------------------
  def plot_lines_and_points(self,states,points,cvec,outfile='maooam-3d',plot_title='3D-slice of MAOOAM attractor',xidx=[0,1,2]):

    nr,nc = np.shape(states)
    x = states[:,xidx[0]]
    y = states[:,xidx[1]]
    z = states[:,xidx[2]]

    nr,nc = np.shape(states)
    xp = points[:,xidx[0]]
    yp = points[:,xidx[1]]
    zp = points[:,xidx[2]]

    trace0 = go.Scatter3d(
      x=x, y=y, z=z,
      mode='lines-and-markers',
      marker=dict(
          size=1,
          color='rgb(0,0,0)',
          opacity=0.5,
      ),
      line=dict(
          color='rgb(0,0,0)',
          width=1
      )
    )

    trace1 = go.Scatter3d(
      x=xp, y=yp, z=zp,
      mode='markers',
      marker=dict(
          size=2,
          color=cvec,
          colorscale='Viridis',
          opacity=0.5,
#         symbol='circle-open',
      ),
#     line=dict(
#         color=cvec,
#         width=1
#     )
    )

    data = [trace0, trace1]

    layout = dict(
      width=1200,
      height=700,
      autosize=False,
      title=plot_title,
      scene=dict(
          xaxis=dict(
              gridcolor='rgb(255, 255, 255)',
              zerolinecolor='rgb(255, 255, 255)',
              showbackground=True,
              backgroundcolor='rgb(230, 230,230)'
          ),
          yaxis=dict(
              gridcolor='rgb(255, 255, 255)',
              zerolinecolor='rgb(255, 255, 255)',
              showbackground=True,
              backgroundcolor='rgb(230, 230,230)'
          ),
          zaxis=dict(
              gridcolor='rgb(255, 255, 255)',
              zerolinecolor='rgb(255, 255, 255)',
              showbackground=True,
              backgroundcolor='rgb(230, 230,230)'
          ),
          camera=dict(
              up=dict(
                  x=0,
                  y=0,
                  z=1
              ),
              eye=dict(
                  x=-1.7428,
                  y=1.0707,
                  z=0.7100,
              )
          ),
          aspectratio = dict( x=1, y=1, z=1 ),
          aspectmode = 'manual'
      ),
    )

    fig = dict(data=data, layout=layout)
    py.plot(fig, filename=outfile, validate=False)

  #-----------------------------------------------------------------------------
  # Plot model output lines and lines
  #-----------------------------------------------------------------------------
  def plot_lines_and_lines(self,states1,states2,cvec,outfile='maooam-3d',plot_title='3D-slice of MAOOAM attractor',name1='trajectory 1', name2='trajectory 2', xidx = [0,1,2]):

    print('states1 = ')
    print(states1)
    print('states2 = ')
    print(states2)
    print('difference = ')
    print(states2 - states1)

    x1 = states1[:,xidx[0]]
    y1 = states1[:,xidx[1]]
    z1 = states1[:,xidx[2]]

    x2 = states2[:,xidx[0]]
    y2 = states2[:,xidx[1]]
    z2 = states2[:,xidx[2]]

    trace1 = go.Scatter3d(
      x=x1, y=y1, z=z1,
      scene='scene1',
      mode='lines',
      line=dict(
          color='rgb(0,0,0)',
          width=2
      ),
      name=name1
    )

    trace2 = go.Scatter3d(
      x=x2, y=y2, z=z2,
      scene='scene2',
      mode='lines',
      line=dict(
          color='rgb(205,12,24)',
          width=2,
      ),
      marker=dict(
          size=4,
          color='rgb(205,12,24)',
      ),
      name=name2
    )

    trace3 = go.Scatter3d(
      x=x2-x1, y=y2-y1, z=z2-z1,
      scene='scene3',
      mode='lines',
      line=dict(
          colorscale='Viridis',
          color=cvec,
          width=2
      ),
      name='difference'
    )

    data = [trace1, trace2]

    eye=np.array([-1.7428,1.0707,0.7100])*1.5
    scene = dict(
          xaxis=dict(
              gridcolor='rgb(255, 255, 255)',
              zerolinecolor='rgb(255, 255, 255)',
              showbackground=True,
              backgroundcolor='rgb(230, 230,230)'
          ),
          yaxis=dict(
              gridcolor='rgb(255, 255, 255)',
              zerolinecolor='rgb(255, 255, 255)',
              showbackground=True,
              backgroundcolor='rgb(230, 230,230)'
          ),
          zaxis=dict(
              gridcolor='rgb(255, 255, 255)',
              zerolinecolor='rgb(255, 255, 255)',
              showbackground=True,
              backgroundcolor='rgb(230, 230,230)'
          ),
          camera=dict(
              up=dict(
                  x=0,
                  y=0,
                  z=1
              ),
              eye=dict(
                  x=eye[0],
                  y=eye[1],
                  z=eye[2]
              )
          ),
          aspectratio = dict( x=1, y=1, z=1 ),
          aspectmode = 'manual'
    )

    layout = dict(
      width=1200,
      height=800,
      autosize=False,
      title=plot_title,
      scene=scene
    )

    fig = tools.make_subplots(rows=1,cols=3,specs=[[{'is_3d':True}, {'is_3d':True},  {'is_3d':True}]])

    # Adding surfaces to subplots
    fig.append_trace(trace1,1,1)
    fig.append_trace(trace2,1,2)
    fig.append_trace(trace3,1,3)

    fig['layout'].update(layout)
    fig['layout']['scene1'].update(scene)
    fig['layout']['scene2'].update(scene)
    fig['layout']['scene3'].update(scene)

##  fig = dict(data=data, layout=layout)
    py.plot(fig, filename=outfile, validate=False, auto_open=False)
