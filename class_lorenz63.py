import numpy as np
from scipy.integrate import odeint
import plotly.offline as py
import plotly.graph_objs as go
from plotly import tools
from copy import deepcopy

#------------------------------------------------------------------
# Define functions
#------------------------------------------------------------------
def f(state, t, sigma, rho, beta):
  x, y, z = state  # unpack the state vector
  return sigma * (y - x), x * (rho - z) - y, x * y - beta * z  # derivatives

def Ja(state, t, sigma, rho, beta):
  x, y, z = state  # unpack the state vector
  J = [[-sigma, sigma,   0   ],
       [ rho-z,    -1,  -x   ],
       [     y,     x,  -beta]]  
  J = np.matrix(J)
  return J

def Jfd(state0,state1,params):
  x0,y0,z0=state0
  x1,y1,z1=state1
  sigma,rho,beta=params

  # If the states are not different, then use machine epsilon as a perturbation
  dist = np.linalg.norm(state1-state0)
  eps = np.finfo(float).eps
  if (dist < eps):
    x1,x2,x3 = eps*np.ones_like(state1)
    print('state0 = ')
    print(state0)
    print('state1 = ')
    print(state1)

  f0 = f([x0,y0,z0],0,sigma,rho,beta)
  fx = f([x1,y0,z0],0,sigma,rho,beta)
  fy = f([x0,y1,z0],0,sigma,rho,beta)
  fz = f([x0,y0,z1],0,sigma,rho,beta)

  # Row 1:
  df1dx = (fx[0] - f0[0])/(x1 - x0)
  df1dy = (fy[0] - f0[0])/(y1 - y0)
  df1dz = (fz[0] - f0[0])/(z1 - z0)

  # Row 2:
  df2dx = (fx[1] - f0[1])/(x1 - x0)
  df2dy = (fy[1] - f0[1])/(y1 - y0)
  df2dz = (fz[1] - f0[1])/(z1 - z0)

  # Row 3:
  df3dx = (fx[2] - f0[2])/(x1 - x0)
  df3dy = (fy[2] - f0[2])/(y1 - y0)
  df3dz = (fz[2] - f0[2])/(z1 - z0)

  J  = np.array([[ df1dx, df1dy, df1dz ],
                 [ df2dx, df2dy, df2dz ],
                 [ df3dx, df3dy, df3dz ]])

  return J


def Jfda(state0,state1,params):
  x0,y0,z0=state0
  x1,y1,z1=state1
  sigma,rho,beta=params

  # If the states are not different, then use machine eps as a perturbation
  dist = np.linalg.norm(state1-state0)
  eps = np.finfo(float).eps
  if (dist < eps):
    x1,x2,x3 = 2.0*eps*np.ones_like(state1)

  f0 = f([x0,y0,z0],0,sigma,rho,beta)
  fx = f([x1,y0,z0],0,sigma,rho,beta)
  fy = f([x0,y1,z0],0,sigma,rho,beta)
  fz = f([x0,y0,z1],0,sigma,rho,beta)

  # Row 1:
  df1dx = -sigma
# df1dx = (fx[0] - f0[0])/(x1 - x0)
  df1dy =  sigma
# df1dy = (fy[0] - f0[0])/(y1 - y0)
  df1dz =  0
# df1dz = (fz[0] - f0[0])/(z1 - z0)

  # Row 2:
  df2dx = rho - (z0+z1)/2.0
# df2dx = (fx[1] - f0[1])/(x1 - x0)
  df2dy = -1
# df2dy = (fy[1] - f0[1])/(y1 - y0)
  df2dz = -(x0+x1)/2.0
# df2dz = (fz[1] - f0[1])/(z1 - z0)

  # Row 3:
  df3dx = (y0+y1)/2.0
# df3dx = (fx[2] - f0[2])/(x1 - x0)
  df3dy = (x0+x1)/2.0
# df3dy = (fy[2] - f0[2])/(y1 - y0)
  df3dz = -beta
# df3dz = (fz[2] - f0[2])/(z1 - z0)

  J  = np.array([[ df1dx, df1dy, df1dz ],
                 [ df2dx, df2dy, df2dz ],
                 [ df3dx, df3dy, df3dz ]])

  return J 

class lorenz63:

  #------------------------------------------------------------------
  # Initialize model parameters
  #------------------------------------------------------------------
  def __init__(self, sigma=10.0, rho=28.0, beta=8.0/3.0):
    self.sigma = sigma
    self.rho = rho
    self.beta = beta
    self.params = [self.sigma,self.rho,self.beta]

  #------------------------------------------------------------------
  # Run model
  #------------------------------------------------------------------
  def run(self, state0, t):
    states = odeint(f, state0, t, args=(self.sigma, self.rho, self.beta))
    return states

  #------------------------------------------------------------------
  # Compute approximate TLM with I + Df(x0)*dt
  #------------------------------------------------------------------
  def compute_TLMa(self, states, t):

    print('states = ')
    print(states)

    print('times = ')
    print(t)

    nr,nc = np.shape(states)
    I = np.identity(nc)

    # Compute Jacobian / linear propagator for each timestep
    sigma,rho,beta = self.params
    maxit = len(t)
    Mhist=[]
    for i in range(maxit):
      if (i < maxit-1):
        dt = t[i+1] - t[i]
      else:
        dt = t[-1] - t[-2]

      # Evaluate Jacobian
      Df = Ja(states[i,:], t[i], sigma, rho, beta)

#     print('Df = ')
#     print(Df)

      # Compute approximate linear propagator
      # (Note this is a poor approximation of the matrix integral, 
      # we would prefer a geometric integrator, e.g. the Magnus expansion)
      M = I + Df*dt  
      Mhist.append(deepcopy(M))

    return Mhist

  #------------------------------------------------------------------
  # Compute approximate Jacobian with finite differences
  #------------------------------------------------------------------
  def compute_Jfd(self, states, t):

    print('states = ')
    print(states)

    print('times = ')
    print(t)

    # Compute Jacobian / linear propagator for each timestep
    Jhist=[]
    for i in range(len(t)-1):
      J = Jfd(states[i,:],states[i+1,:],self.params)
      print('J = ')
      print(J)
      Jhist.append(J)

    return Jhist

  #------------------------------------------------------------------
  # Compute approximate Jacobian with partial analytic and partial finite differences
  #------------------------------------------------------------------
  def compute_Jfda(self, states, t):

    print('states = ')
    print(states)

    # Compute TLM for each timestep
    Jhist=[]
    for i in range(len(t)-1):
      J = Jfda(states[i,:],states[i+1,:],self.params)
      Jhist.append(J)

    return Jhist

  #------------------------------------------------------------------
  # Plot model output
  #------------------------------------------------------------------
  def plot(self, states,cvec,outfile='l63-3d',plot_title='Lorenz 63 attractor'):

    nr,nc = np.shape(states)
    x = states[:,0]
    y = states[:,1]
    z = states[:,2]

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

  #------------------------------------------------------------------
  # Plot model output lines and points
  #------------------------------------------------------------------
  def plot_lines_and_points(self,states,points,cvec,outfile='l63-3d',plot_title='Lorenz 63 attractor'):

    nr,nc = np.shape(states)
    x = states[:,0]
    y = states[:,1]
    z = states[:,2]

    nr,nc = np.shape(states)
    xp = points[:,0]
    yp = points[:,1]
    zp = points[:,2]

    trace0 = go.Scatter3d(
      x=x, y=y, z=z,
      mode='lines-and-markers',
      marker=dict(
          size=1,
          color=cvec,
          colorscale='Viridis',
          opacity=0.5,
      ),
      line=dict(
          color=cvec,
          width=1
      )
    )

    trace1 = go.Scatter3d(
      x=xp, y=yp, z=zp,
      mode='markers',
      marker=dict(
          size=3,
          color=cvec,
          colorscale='Viridis',
          symbol='circle-open',
      ),
#     line=dict(
#         color=cvec,
#         width=1
#     )
    )

    data = [trace0, trace1]

    layout = dict(
      width=800,
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

  #------------------------------------------------------------------
  # Plot model output lines and lines
  #------------------------------------------------------------------
  def plot_lines_and_lines(self,states1,states2,cvec,outfile='l63-3d',plot_title='Lorenz 63 attractor',name1='trajectory 1', name2='trajectory 2'):

    print('states1 = ')
    print(states1)
    print('states2 = ')
    print(states2)
    print('difference = ')
    print(states2 - states1)

    x1 = states1[:,0]
    y1 = states1[:,1]
    z1 = states1[:,2]

    x2 = states2[:,0]
    y2 = states2[:,1]
    z2 = states2[:,2]

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
    py.plot(fig, filename=outfile, validate=False)


