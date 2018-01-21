import numpy as np
import pickle

class state_vector:
  def __init__(self,params=[0],x0=[0],t=[0],name='uninitialized'):
    self.tdim = np.size(t)
    self.xdim = np.size(x0)
    self.params = params
    self.x0 = x0
    self.t = t
    self.name = name
    self.trajectory = np.zeros((self.tdim,self.xdim))
    self.clim_mean = [0,0,0]
    self.clim_std  = [0,0,0]
    self.Jhist = []
    self.Mhist = []
    self.Qhist = []
    self.Rhist = []
    self.hist_ainc = 0
    self.hist_dtau = 0
    self.hist_idx = []
    self.rescale_interval = 1
    self.TLM = []

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
    print('Climatological Mean:')
    print(self.clim_mean)
    print('Climatological Standard Deviation:')
    print(self.clim_std)
    return self.name

  def setName(self,name):
    self.name = name

  def getTrajectory(self):
    return self.trajectory

  def setTrajectory(self,states):
    self.trajectory = states

    # Compute the climatological mean and variance
    x_avg = np.zeros(self.xdim)
    x_std = np.zeros(self.xdim)
    for i in range(self.xdim):
      x_avg[i] = np.mean(states[:,i])
      x_std[i] = np.std(states[:,i])
    self.clim_mean = x_avg
    self.clim_std = x_std

  def getClimMean(self):
    return self.clim_mean

  def getClimStd(self):
    return self.clim_std

  def getTLM(self):
    return self.TLM

  def setTLM(self,TLM):
    self.TLM = TLM

  def getJhist(self):
    return self.Jhist

  def setJhist(self,Jhist):
    self.Jhist = Jhist

  def getMhist(self):
    return self.Mhist

  def setMhist(self,Mhist):
    self.Mhist = Mhist

  def getM2hist(self):
    return self.M2hist

  def setM2hist(self,M2hist):
    self.M2hist = M2hist

  def getQhist(self):
    return self.Qhist

  def setQhist(self,Qhist):
    self.Qhist = Qhist

  def getRhist(self):
    return self.Rhist

  def setRhist(self,Rhist):
    self.Rhist = Rhist

  def getLEs(self):
    return self.LEs

  def setLEs(self,LEs):
    self.LEs = LEs

  def getTimes(self):
    return self.t

  def save(self,outfile):
    with open(outfile,'wb') as output:
      pickle.dump(self,output,pickle.HIGHEST_PROTOCOL)

  def load(self,infile):
    with open(infile,'rb') as input:
      sv = pickle.load(input)
    return sv
