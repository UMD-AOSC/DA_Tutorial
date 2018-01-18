import numpy as np
import pickle

class obs_data:
  def __init__(self,t=[0],pos=[0],val=[0],err=[0],bias=[],xt=[0],name='uninitialized'):
    tdim = 0
    xdim = 0
    odim = 0
    self.name = name
    self.t = t
    self.pos = np.array(pos)
    self.val = np.array(val)
    self.err = np.array(err)
    self.bias = np.array(bias)
    self.hx = np.array(val)
    self.xt = np.array(xt)
    self.dep = []
    self.climvar = []

  def __str__(self):
    print(self.name)
    print('Obs Error:')
    print(self.err)
    print('Positions:')
    print(self.pos)
    print('Time interval:')
    print(self.t)
    print('Values:')
    print(self.val)
    return self.name

  def getVal(self):
    return self.val

  def getErr(self):
    return self.err

  def getPos(self):
    return self.pos

  def getDep(self):
    return self.dep

  def getHx(self):
    return self.hx

  def setVal(self,val):
    self.val = np.array(val)

  def setErr(self,err):
    self.err = np.array(err)

  def setPos(self,pos):
    self.pos = np.array(pos)

  def setDep(self,dep):
    self.dep = np.array(dep)

  def setHx(self,hx):
    self.hx = np.array(hx)

  def save(self,outfile):
    with open(outfile,'wb') as output:
      pickle.dump(self,output,pickle.HIGHEST_PROTOCOL)

  def load(self,infile):
    with open(infile,'rb') as input:
      sv = pickle.load(input)
    return sv
