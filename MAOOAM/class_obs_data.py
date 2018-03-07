import numpy as np
import pickle

class obs_data:
  def __init__(self,t=[0],pos=[0],val=[0],err=[0],bias=[0],xt=[0],name='uninitialized',mu_init=[],sigma_init=[]):
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
    self.mu_init = mu_init
    self.sigma_init = sigma_init

  def __str__(self):
    print(self.name)
    print('mu_init = ', self.mu_init)
    print('sigma_init = ', self.sigma_init)
    print('Obs Error:')
    print(self.err)
    print('Positions:')
    print(self.pos)
#   print('Time interval:')
#   print(self.t)
    print('Values:')
    print(self.val)
    return self.name

  def getVal(self):
    return np.asarray(self.val)

  def getErr(self):
    return self.err

  def getPos(self):
    return self.pos.astype(int)

  def getDep(self):
    return self.dep

  def getHx(self):
    return self.hx

  def setVal(self,val):
    self.val = np.array(val)

  def setErr(self,err):
    self.err = np.array(err)

  def setPos(self,pos):
    self.pos = np.array(pos).astype(int)

  def setDep(self,dep):
    self.dep = np.array(dep)

  def setHx(self,hx):
    self.hx = np.array(hx)

  def reduceDim(self,dims_to_keep):
    if (len(dims_to_keep)==1):
      self.val = self.val[:,dims_to_keep[0]:dims_to_keep[0]+1]
      self.err = self.err[:,dims_to_keep[0]:dims_to_keep[0]+1]
      self.pos = self.pos[:,dims_to_keep[0]:dims_to_keep[0]+1]
      self.dep = self.dep[:,dims_to_keep[0]:dims_to_keep[0]+1]
      self.hx = self.hx[:,dims_to_keep[0]:dims_to_keep[0]+1]
    else:
      self.val = self.val[:,dims_to_keep]
      self.err = self.err[:,dims_to_keep]
      self.pos = self.pos[:,dims_to_keep]
      self.dep = self.dep[:,dims_to_keep]
      self.hx = self.hx[:,dims_to_keep]
  
  def fillDim(self,dims_to_fill,fillValue):
    nr = self.tdim
    nc = self.xdim
    val = np.ones((nr,nc))*fillValue
    err = np.ones((nr,nc))*fillValue
    dep = np.ones((nr,nc))*fillValue
    hx  = np.ones((nr,nc))*fillValue

    val[:,self.pos] = self.val
    err[:,self.pos] = self.err
    dep[:,self.pos] = self.dep
    hx[:,self.pos]  = self.hx

    self.val = val
    self.err = err
    self.dep = dep
    self.hx = hx

    for i in range(self.tdim): 
      self.pos[i,:] = range(self.xdim)
    
  def save(self,outfile):
    with open(outfile,'wb') as output:
      pickle.dump(self,output,pickle.HIGHEST_PROTOCOL)

  def load(self,infile):
    with open(infile,'rb') as input:
      sv = pickle.load(input)
    return sv
