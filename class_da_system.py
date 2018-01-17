import numpy as np
import scipy as sp
from class_state_vector import state_vector
from class_obs_data import obs_data

class da_system:

  def __init__(self,x0=[0],yo=[0],t0=0,dt=0):
    self.xdim = np.size(x0)
    self.ydim = np.size(yo)
    self.x0 = x0
    self.t0 = t0
    self.dt = dt
    self.X0 = x0
    self.t = t0
    self.B = np.matrix(np.identity(self.xdim))
    self.R = np.matrix(np.identity(self.ydim))
    self.H = np.matrix(np.identity(self.xdim))
    self.Ht = (self.H).transpose()
#   self.SqrtB = 

  def __str__(self):
    print('xdim = ', self.xdim)
    print('ydim = ', self.ydim)
    print('x0 = ', self.x0)
    print('t0 = ', self.t0)
    print('dt = ', self.dt)
    print('t  = ', self.t)
    print('B = ')
    print(self.B)
    print('R = ')
    print(self.R)
    return 'type::da_system'

  def setMethod(self,method):
    self.method = method

  def update(self,B=[0],R=[0],H=[0],t=[0],x0=[0]):
    # Update the state of the DA system for a new cycle
    self.B = B
    self.R = R
    self.H = H
    self.Ht = H.Transpose()
    self.t = t
    self.x0 = x0

  def getB(self):
    return self.B

  def setB(self,B):
    self.B = np.matrix(B)
    nr,nc = np.shape(B)
    self.xdim = nr

  def getR(self):
    return self.R

  def setR(self,R):
    self.R = np.matrix(R)
    nr,nc = np.shape(R)
    self.ydim = nr

  def getH(self):
    return self.H

  def setH(self,H):
    self.H = np.matrix(H)
    self.Ht = np.transpose(self.H)
    nr,nc = np.shape(H)

    # Verify that H is ydim x xdim
    if (nr != self.ydim):
      error('H must be ydim x xdim, but instead H is %d x %d'%(nr,nc))
    if (nc != self.xdim):
      error('H must be ydim x xdim, but instead H is %d x %d'%(nr,nc))

  def compute_analysis(self,xb,yo,params=[0]):
    method = self.method
    if method == 'skip':
      xa = xb 
    elif method == 'nudging':
      xa = self.nudging(xb,yo)
    elif method == 'OI':
      xa = self.OI(xb,yo)
    elif method == '3DVar' or method == '3D-Var':
      xa = self._3DVar(xb,yo) 
#   elif method == '4DVar' or method == '4D-Var':
#     xa = self._4DVar(xb,yo)
#   elif method == 'PF':
#     xa = self.PF(xb,yo)
#   elif method == 'ETKF' or method == 'EnKF':
#     xa = self.ETKF(xb,yo)
#   elif method == '4DEnVar':
#     xa = self._4DEnVar(xb,yo)
#   elif method == '4DETKF':
#     xa = self._4DETKF(xb,yo) 
    else:
      print('compute_analysis:: Unrecognized DA method.')
    return xa

# def init_ens(self,sigma=0.1,nens=2):
#   xdim = size(self.x0)
#   self.X0 = np.matlib.repmat(self.x0, 1, nens) + np.random.normal(mu,sigma,(xdim,nens))

# def init_4D(self):
#   xdim = size(self.x0)

# def init_4DEns(self):
#   xdim = size(self.x0)

#---------------------------------------------------------------------------------------------------
  def nudging(self,xb,yo):
#---------------------------------------------------------------------------------------------------
    # Use observations at predefined points to drive the model system to the observed nature system
    const = np.diagonal(self.B)
    xa = xb + const*(yo - xb) 

    C = np.diag(const)
    Hl = self.H
    Ht = self.Ht
    KH = Ht*C*Hl

    return xa,KH
        
#---------------------------------------------------------------------------------------------------
  def OI(self,xb,yo):
#---------------------------------------------------------------------------------------------------
    xb = np.matrix(xb).flatten().T
    yo = np.matrix(yo).flatten().T

    # Use explicit expression to solve for the analysis
    print(self)
    H = self.H
    Ht = self.Ht
    B = self.B
    R = self.R

    KH = B*Ht*np.linalg.inv(H*B*Ht + R)*H

    print('KH = ')
    print(KH)

    Hxb = H*xb.flatten().T
    xa = xb + KH*(yo - Hxb)

    return xa.A1, KH


#---------------------------------------------------------------------------------------------------
  def _3DVar(self,xb,yo):
#---------------------------------------------------------------------------------------------------
    # Use minimization algorithm to solve for the analysis
    xdim = self.xdim
    Hl = self.H
    Ht = self.Ht
    B = self.B
    R = self.R
    Rinv = np.linalg.inv(R)

    # make inputs column vectors
    xb = np.matrix(xb).flatten().T
    yo = np.matrix(yo).flatten().T

    # 'preconditioning with B'
    I = np.identity(xdim)
    BHt = B*Ht
    BHtRinv = BHt*Rinv
    A = I + BHtRinv*Hl
    b1 = xb + BHtRinv*yo

    # Use minimization algorithm to minimize cost function:
    xa,ierr = sp.sparse.linalg.cg(A,b1,x0=xb,tol=1e-05,maxiter=1000) 
#   xa,ierr = sp.sparse.linalg.bicgstab(A,b1,x0=np.zeros_like(b1),tol=1e-05,maxiter=1000)
#   try: gmres, 
#   xa = sp.optimize.minimize(fun, x0, args=(), method=None, jac=None, hess=None, hessp=None, bounds=None, constraints=(), tol=None, callback=None, options=None)

    # Compute KH:
    HBHtPlusR_inv = np.linalg.inv(Hl*BHt + R)
    KH = BHt*HBHtPlusR_inv*Hl

    return xa, KH


#---------------------------------------------------------------------------------------------------
# def ETKF(self,Xb,yo):
#---------------------------------------------------------------------------------------------------
#   # Use ensemble of states to estimate background error covariance
    xdim = self.xdim
    Yf = self.H*Xb
    R = self.R
    Yb = np.zeros([nobs,kdim])
    for i in range(kdim):
      Yb[:,i] = H*Xb[:,i]
    
    # Convert ensemble members to perturbations
    xm = np.mean(Xb,axis=1)
    ym = np.mean(Yb,axis=1)
    Xb = Xb - xm[:,np.newaxis]
    Yb = Yb - ym[:,np.newaxis]

    # Compute R^{-1}
    Rinv = np.linalg.inv(R)

    # Compute the weights
    Ybt = np.transpose(Yb)
    C = np.dot(Ybt,Rinv)

    I = np.identity(kdim)
    rho = 1.0
    eigArg = (kdim-1)*I/rho + np.dot(C,Yb)
    lamda,P = np.linalg.eigh(eigArg)
    Linv = np.diag(1.0/lamda)
    PLinv = np.dor(P,Linv)
    Pt = np.transpose(P)
    Pa = np.dot(PLinv,Pt)

    Linvsqrt = np.diag(1/np.sqrt(lamda))
    PLinvsqrt = np.dor(P,Linvsqrt)
    Wa = np.sqrt(kdim-1) * np.dot(PLinvsqrt,Pt)

    d = yo - ym
    Cd = np.dot(C,d)
    wm = np.dot(Pa,Cd)

    Wa = Wa + wm[:,np.newaxis]

    # Add the same mean vector wm to each column
    Xa = np.dot(Xb,Wa) + xm[:,np.newaxis]

    # Compute KH:
    Hl = self.H
    Pb = np.dot(Xb,np.transpose(Xb))
    Ht = np.transpose(H)
    PbHt = np.dot*(Pb,Ht)
    HPbHtPlusR_inv = np.linalg.inv(np.dot(Hl,PbHt) + R)
    K = np.dot(PbHt,HPbHtPlusR_inv)
    KH = np.dot(K,Hl)
    
    
    return Xa, KH


#---------------------------------------------------------------------------------------------------
# def _4DVar(self,xb_4d,yo_4d):
#---------------------------------------------------------------------------------------------------
#   # Use minimization over a time window to solve for the analysis
    
    B = self.B
    R = self.R
    xdim = self.xdim
    Xb = np.matrix(np.atleast_2d(Xb))
    Yo = np.matrix(np.atleast_2d(Yo))
    Yp = np.matrix(np.atleast_2d(Yp))



#---------------------------------------------------------------------------------------------------
# def 4DEnVar(self,Xb_4d,yo_4d):
#---------------------------------------------------------------------------------------------------
#   # Use ensemble of states over a time window to estimate temporal correlations
#   xdim = size(self.x0)


#---------------------------------------------------------------------------------------------------
# def 4DETKF(self,Xb_4d,yo_4d):
#---------------------------------------------------------------------------------------------------
#   # Use ensemble of states over a time window to estimate temporal correlations
#   xdim = size(self.x0)


#---------------------------------------------------------------------------------------------------
# def PF(self,Xb,yo):
#---------------------------------------------------------------------------------------------------
#   # Use an ensemble of states to estimate a multidimensional probability distribution
#   xdim = size(self.x0)


#---------------------------------------------------------------------------------------------------
# def Hybrid(self,Xb,yo):
#---------------------------------------------------------------------------------------------------
#   # Use a hybrid method to compute the analysis
#   xdim = size(self.x0)
