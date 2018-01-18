import numpy as np
import pickle

class state_vector:
  def __init__(self,params=[0],x0=[0],t=[0],name='uninitialized'):
    tdim = np.size(t)
    xdim = np.size(x0)
    self.params = params
    self.x0 = x0
    self.t = t
    self.name = name
    self.trajectory = np.zeros((tdim,xdim))

  def __str__(self):
    print(self.name)
    print('Parameters:')
    print(self.params)
    print('Initial conditions:')
    print(self.x0)
    print('Time array:')
    print(self.t)
    print('Trajectory:')
    print(self.trajectory)
    return self.name

  def getTrajectory(self):
    return self.trajectory

  def setTrajectory(self,states):
    self.trajectory = states

  def getTLM(self):
    return self.TLM

  def setTLM(self,TLM):
    self.TLM = TLM

  def getTimes(self):
    return self.t

  def save(self,outfile):
    with open(outfile,'wb') as output:
      pickle.dump(self,output,pickle.HIGHEST_PROTOCOL)

  def load(self,infile):
    with open(infile,'rb') as input:
      sv = pickle.load(input)
    return sv
