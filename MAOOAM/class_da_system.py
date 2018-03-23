import numpy as np
import scipy as sp
from class_state_vector import state_vector
from class_obs_data import obs_data
import numpy.matlib
import pickle
from copy import deepcopy
from module_obs_network import get_h_full_coverage

#===============================================================================
class da_system:
#===============================================================================

  #-----------------------------------------------------------------------------
  def __init__(self,x0=[],yo=[],t0=0,dt=0,alpha=0.5,state_vector=[],obs_data=[],acyc_step=10):
  #-----------------------------------------------------------------------------
    self.xdim = np.size(x0)
    self.ydim = np.size(yo)
    self.edim = 1
    self.x0 = x0
    self.t0 = t0
    self.dt = dt
    self.X0 = x0
    self.t = t0
    self.acyc_step = acyc_step
    self.dtau = dt*acyc_step
    self.fcst_step = acyc_step
    self.fcst_dt = dt
    self.maxit = 0
    self.B = np.matrix(np.identity(self.xdim))
    self.R = np.matrix(np.identity(self.ydim))
    # self.H = np.matrix(get_h_full_coverage())
    self.H = np.matrix(np.identity(self.xdim))
    self.Ht = (self.H).transpose()
    self.alpha = alpha
    self.SqrtB = []
    self.state_vector = state_vector
    self.obs_data = obs_data
    self.method = ''
    self.KH = []
    self.khidx = []
    self.edim = 0
    self.das_bias_init = 0
    self.das_sigma_init = 1
    self.outer_loops = 1

  #-----------------------------------------------------------------------------
  def __str__(self):
  #-----------------------------------------------------------------------------
    print('xdim = ', self.xdim)
    print('ydim = ', self.ydim)
    print('x0 = ', self.x0)
    print('t0 = ', self.t0)
    print('dt = ', self.dt)
    print('t  = ', self.t)
    print('acyc_step = ', self.acyc_step)
    print('dtau = ', self.dtau)
    print('fcst_step = ', self.fcst_step)
    print('fcst_dt  = ', self.fcst_dt)
    print('B = ')
    print(self.B)
    print('R = ')
    print(self.R)
    print('H = ')
    print(self.H)
    print('state_vector = ')
    print(self.state_vector)
    print('obs_data = ')
    print(obs_data)
    print('method = ')
    print(self.method)
    return 'type::da_system'

  #-----------------------------------------------------------------------------
  # Set / Get methods
  #-----------------------------------------------------------------------------

  def setMethod(self,method):
    self.method = method

  def getMethod(self):
    return self.method

  def setStateVector(self,sv):
    self.state_vector = sv

  def setObsData(self,obs):
    self.obs_data = obs

  def getStateVector(self):
    return self.state_vector

  def getObsData(self):
    return self.obs_data

  def getC(self):
    return self.C

  def setC(self,C):
    self.C = np.matrix(C)

  def getB(self):
    return self.B

  def setB(self,B):
    self.B = np.matrix(B)
    nr,nc = np.shape(B)
    self.xdim = nr
    self.SqrtB = sp.linalg.sqrtm(self.B)

  def setSqrtB(self,X):
    self.SqrtB = np.matrix(X)
    self.B = self.SqrtB*self.SqrtB.T
    nr,nc = np.shape(X)
    self.xdim = nr
    self.edim = nc

  def getR(self):
    return self.R

  def setR(self,R):
    self.R = np.matrix(R)
#   nr,nc = np.shape(R)
#   self.ydim = nr
    self.Rinv = np.linalg.inv(R)

  def getH(self):
    return self.H

  def setH(self,H):
    self.H = np.matrix(H)
    self.Ht = np.transpose(self.H)
#   nr,nc = np.shape(H)

    # Verify that H is ydim x xdim
#   if (nr != self.ydim):
#     error('H must be ydim x xdim, but instead H is %d x %d'%(nr,nc))
#   if (nc != self.xdim):
#     error('H must be ydim x xdim, but instead H is %d x %d'%(nr,nc))

  def getKH(self):
    return self.KH, self.khidx

  def setKH(self,KH,khidx):
    self.KH = KH
    self.khidx = khidx

  #-----------------------------------------------------------------------------
  # Support functions
  #-----------------------------------------------------------------------------

  def reduceYdim(self,yp):
#   print('reduceYdim:')
#   print('yp = ', yp)
    self.ydim = len(yp)
    self.setH(self.H[yp,:])
    R = self.R
    R = R[yp,:]
    R = R[:,yp]
    self.setR(R)

  def expandToList(self, M, n):
    # Copy one matrix 'n' times into a list of length 'n'
    Mlist = []
    for i in range(n):
      Mlist.append(deepcopy(M))
    return Mlist    

  def convertList2Stacked(self,Xlist):
    # Setup the augmented ensemble of state vectors by concatenating all timesteps
    nl = len(Xlist)
    nr,nc = np.shape(Xlist[0])
    X = np.zeros((nr*nl,nc))
    for i in range(nl):
      X[i*nr:(i+1)*nr,:] = Xlist[i]
    return X

  def expandStackedMatrix(self,Alist):
    # Takes a list of matrices (e.g. error covariance matrices) and builds a
    # single large matrix
    verbose = True
    nl = len(Alist)
    nr,nc = np.shape(Alist[0])
    print('nr = ', nr)
    print('nc = ', nc)
    A = np.zeros((nr*nl,nc*nl))
    for i in range(nl): 
      A[i*nr:(i+1)*nr,i*nc:(i+1)*nc] = Alist[i]
    return A

  # (Unused)
  def update(self,B=[0],R=[0],H=[0],t=[0],x0=[0]):
    # Update the state of the DA system for a new cycle
    self.B = B
    self.R = R
    self.H = H
    self.Ht = H.Transpose()
    self.t = t
    self.x0 = x0

  def initEns(self,x0,mu=0,sigma=0.1,edim=4):
    # Initialize a set of ensemble members
    xdim = len(x0)
    x0 = np.matrix(x0).flatten().T
    mu = np.matrix(mu).flatten().T
    Xrand = np.random.normal(0,sigma,(xdim,edim))
    Xrand = np.matrix(Xrand)
#   print('Xrand = ')
#   print(Xrand)
    # Remove mean to make sure it is properly centered at 0
    # (add bias if included)
    rand_mean = np.mean(Xrand,axis=1) - mu
#   print('rand_mean = ')
#   print(rand_mean)
    rmat = np.matlib.repmat(rand_mean,1,edim)
    Xrand = Xrand - rmat
#   print('xrand = ')
#   print(xrand)
    # add perturbations to x0
    rmat = np.matlib.repmat(x0, 1, edim)
#   print('rmat = ')
#   print(rmat)
    X0 = np.matrix(Xrand + rmat)
    return X0

# def init4DEns(self):
#   xdim = size(self.x0)

  #-----------------------------------------------------------------------------
  def compute_analysis(self,xb,yo,params=[0]):
  #-----------------------------------------------------------------------------
    # (params needed for 4D-Var)
#   method = self.method
    method = str.lower(self.method)
    if method == str.lower('skip'):
      xa,KH = xb,np.zeros(len(xb))
    elif method == str.lower('nudging'):
      xa,KH = self.nudging(xb,yo)
    elif method == str.lower('OI'):
      xa,KH = self.OI(xb,yo)
    elif method == str.lower('3DVar') or method == str.lower('3D-Var'):
      xa,KH = self._3DVar(xb,yo) 
    elif method == str.lower('ETKF') or method == str.lower('EnKF'):
      xa,KH = self.ETKF(xb,yo)
    elif method == str.lower('PF'):
      xa,KH = self.PF(xb,yo)
    elif method == str.lower('Hybrid'):
      xa,KH = self.HybridGain(xb,yo) 
    elif method == str.lower('4DVar') or method == str.lower('4D-Var'):
      xa,KH = self._4DVar_outerLoop(xb,yo,params)
    elif method == str.lower('4DEnVar'):
      xa,KH = self._4DEnVar(xb,yo,params)
    elif method == str.lower('4DETKF'):
      xa,KH = self._4DETKF(xb,yo,params) 
    else:
      print('compute_analysis:: Unrecognized DA method.')
      print('method = ', method)
      raise SystemExit
    return xa,KH

#---------------------------------------------------------------------------------------------------
# DA Methods:
#---------------------------------------------------------------------------------------------------

  #-------------------------------------------------------------------------------------------------
  def nudging(self,xb,yo):
  #-------------------------------------------------------------------------------------------------
  # Use observations at predefined points to drive the model system to the observed nature system

    verbose = False

    xb = np.matrix(xb).flatten().T
    yo = np.matrix(yo).flatten().T
    Hl = np.matrix(self.H)
    Ht = np.matrix(self.Ht)

    C = self.C
    xa = xb + np.dot(C,Ht)*(yo - Hl*xb) 

    if verbose:
      print('xb = ')
      print(xb)
      print('C = ')
      print(C)
      print('yo = ')
      print(yo)
      print('Hl = ')
      print(Hl)
      print('xa = ')
      print(xa)
      print('Ht = ')
      print(Ht)

    KH = np.dot(np.dot(Hl,C),Ht)

    return xa.A1,KH

        
  #-------------------------------------------------------------------------------------------------
  def OI(self,xb,yo):
  #-------------------------------------------------------------------------------------------------

    verbose = False

    xb = np.matrix(xb).flatten().T
    yo = np.matrix(yo).flatten().T

    # Use explicit expression to solve for the analysis
#   print(self)
    Hl = self.H
    Ht = self.Ht
    B = self.B
    R = self.R

    if verbose:
      print('H = ')
      print(Hl)
      print('R = ')
      print(R)

    # Should be dimensioned: xdim * xdim
    K = B*Ht*np.linalg.inv(Hl*B*Ht + R)
    KH = K*Hl

    if verbose:
      print('KH = ')
      print(KH)

    Hxb = Hl*xb
    xa = xb + K*(yo - Hxb)

    return xa.A1, KH


  #-------------------------------------------------------------------------------------------------
  def _3DVar(self,xb,yo):
  #-------------------------------------------------------------------------------------------------
  # Use minimization algorithm to solve for the analysis

    # make inputs column vectors
    xb = np.matrix(xb).flatten().T
    yo = np.matrix(yo).flatten().T

    # Set parameters
    xdim = self.xdim
    Hl = self.H
    Ht = self.Ht
    B = self.B
    R = self.R
    Rinv = np.linalg.inv(R)

    # 'preconditioning with B'
    I = np.identity(xdim)
    BHt = np.dot(B,Ht)
    BHtRinv = np.dot(BHt,Rinv)
    A = I + np.dot(BHtRinv,Hl)
    b1 = xb + np.dot(BHtRinv,yo)

    # Use minimization algorithm to minimize cost function:
    xa,ierr = sp.sparse.linalg.cg(A,b1,x0=xb,tol=1e-05,maxiter=1000) 
#   xa,ierr = sp.sparse.linalg.bicgstab(A,b1,x0=np.zeros_like(b1),tol=1e-05,maxiter=1000)
#   try: gmres, 

    # Compute KH:
    HBHtPlusR_inv = np.linalg.inv(Hl*BHt + R)
    KH = BHt*HBHtPlusR_inv*Hl

    return xa, KH


  #-------------------------------------------------------------------------------------------------
  def ETKF(self,Xb,yo):
  #-------------------------------------------------------------------------------------------------
  # Use ensemble of states to estimate background error covariance
    verbose=False

    # Make sure inputs are matrix and column vector types
    Xb = np.matrix(Xb)
    yo = np.matrix(yo).flatten().T

    # Get system dimensions
    nr,nc = np.shape(Xb)
    xdim = nr  # == self.xdim
    edim = nc
    ydim = len(yo)

    # Apply observation operator to forecast ensemble
    Hl = self.H
    Yb = np.matrix(np.zeros([ydim,edim]))
    for i in range(edim):
      Yb[:,i] = np.dot(Hl,Xb[:,i])
    
    # Convert ensemble members to perturbations
    xm = np.mean(Xb,axis=1)
    ym = np.mean(Yb,axis=1)
    Xb = Xb - np.matlib.repmat(xm, 1, edim)
    Yb = Yb - np.matlib.repmat(ym, 1, edim)

    # Compute R^{-1}
    R = self.R
    Rinv = np.linalg.inv(R)

    # Compute the weights

    #----
    # stage(4) Now do computations for lxk Yb matrix
    # Compute Yb^T*R^(-1)
    #----
    Ybt = np.transpose(Yb)
    C = np.dot(Ybt,Rinv)

    if verbose:
      print ('C = ')
      print (C)

    #----
    # stage(5) Compute eigenvalue decomposition for Pa
    # Pa = [(k-1)I/rho + C*Yb]^(-1)
    #----
    I = np.identity(edim)
    rho = 1.05 #1.0
    eigArg = (edim-1)*I/rho + np.dot(C,Yb)
    print(np.linalg.eig(Rinv)[0]); import sys; sys.exit(1)

    lamda,P = np.linalg.eigh(eigArg)

    if verbose:
      print ('lamda = ')
      print (lamda)
      print ('P = ')
      print (P)

    Linv = np.diag(1.0/lamda)

    if verbose:
      print ('Linv = ')
      print (Linv)

    PLinv = np.dot(P,Linv)

    if verbose:
      print ('PLinv = ')
      print (PLinv)

    Pt = np.transpose(P)

    if verbose:
      print ('Pt = ')
      print (Pt)

    Pa = np.dot(PLinv, Pt)

    if verbose:
      print ('Pa = ')
      print (Pa) 

    #----
    # stage(6) Compute matrix square root
    # Wa = [(k-1)Pa]1/2
    #----
    Linvsqrt = np.diag(1/np.sqrt(lamda))

    if verbose:
      print ('Linvsqrt = ')
      print (Linvsqrt)

    PLinvsqrt = np.dot(P,Linvsqrt)

    if verbose:
      print ('PLinvsqrt = ')
      print (PLinvsqrt)

    Wa = np.sqrt((edim-1)) * np.dot(PLinvsqrt,Pt)

    if verbose:
      print ('Wa = ')
      print (Wa)

    #----
    # stage(7) Transform back
    # Compute the mean update
    # Compute wabar = Pa*C*(yo-ybbar) and add it to each column of Wa
    #----
    d = yo-ym
    Cd = np.dot(C,d)

    if verbose:
      print ('Cd = ')
      print (Cd)

    wm = np.dot(Pa,Cd) #k x 1

    if verbose:
      print ('wm = ')
      print (wm)

    # Add the same mean vector wm to each column
    Wa = Wa + np.matlib.repmat(wm, 1, edim)

    if verbose:
      print ('Wa = ')
      print (Wa)

    #----
    # stage(8)
    # Compute the perturbation update
    # Multiply Xb (perturbations) by each wa(i) and add xbbar
    #----

    # Add the same mean vector wm to each column
    Xa = np.dot(Xb,Wa) + np.matlib.repmat(xm, 1, edim)

    if verbose:
      print ('Xa = ')
      print (Xa)

    # Compute KH:
    RinvYb = np.dot(Rinv,Yb)
    IpYbtRinvYb = ((edim-1)/rho)*I + np.dot(Ybt,RinvYb)
    IpYbtRinvYb_inv = np.linalg.inv(IpYbtRinvYb)
    YbtRinv = np.dot(Ybt,Rinv)
    K = np.dot( Xb, np.dot(IpYbtRinvYb_inv,YbtRinv) )
    KH = np.dot(K,Hl)
    
    return Xa, KH


  #-------------------------------------------------------------------------------------------------
  def _4DVar_outerLoop(self,xb_4d,yo_4d,params):
  #-------------------------------------------------------------------------------------------------
  # Use minimization over a time window to solve for the analysis
#   
#   xb_4d = np.matrix(np.atleast_2d(xb_4d))
#   yo_4d = np.matrix(np.atleast_2d(yo_4d))
#   nr,nc = np.shape(xb_4d)     ! columns are state vectors at consecutive timesteps
#   xdim = nr
#   tdim = nc
#   B = self.B
#   R_4d = self.R               ! may need list of R matrices if observations are non-uniform over time
#   H_4d = self.H               ! may need list of H matrices if observations are non-uniform over time
#   M_4d = self.M               ! input list of TLMs for each timestep
#
# (work in progress)

    # Using 'incremental form' of Courtier et al. (1994)
    # 'quadratic version' (13), p. 1379

    verbose = False

    #-----------------------------------------------------------------------------
    # 4D-Var Outer loop (runs 1 iteration of the outer loop)
    #-----------------------------------------------------------------------------

    if verbose:
      print ('------------------------- In das._4dvarOuterLoop -----------------------------')

    # Using initial conditions of this cycle's forecast
    x0_4dvar = np.matrix(xb_4d[0,:]).flatten().T

    if verbose:
      print('Initial:')
      print('x0_4dvar = ')
      print(x0_4dvar)

    if verbose:
      print('yo_4d = ')
      print(yo_4d)
  
    if verbose:
      print('xb_4d = ')
      print(xb_4d)
  
    # Get new initial conditions x0 via 4dvar inner loop:
    xb_4d = np.matrix(xb_4d).T
    yo_4d = np.matrix(yo_4d).T
    x0_4dvar,KH = self._4DVar_innerLoop(xb_4d,yo_4d,params)

    if verbose:
      print('Initial:')
      print('x0_4dvar = ')
      print(x0_4dvar)

    #-----------------------------------------------------------
    # Return to rerun the model forecast, and repeat this loop
    #-----------------------------------------------------------
    return x0_4dvar,KH


#   # Collect the final analysis:
#   xa_4dvar = xb_4dvar

#   if verbose:
#     print ('------------------------------------------------------')
#     print ('Forecast :: ', x_forecast)
#     print ('4Dvar analysis :: ', xa_4dvar)
#     print ('xn = ')
#     print (xnat[i,:])
#     print ('Forecast error (absolute) :: ')
#     print (abs(x_forecast - xnat[i,:]))
#     print ('Analysis error (absolute) :: ')
#     print (abs(xa_4dvar - xnat[i,:]))
#     print ('------------------------------------------------------')

#     print('Post-4DVar: stopping on purpose.')
#     exit()
   

    return

  #-------------------------------------------------------------------------------------------------
  def _4DVar_innerLoop(self,xb_4d,yo_4d,params):
  #-------------------------------------------------------------------------------------------------

    verbose = False

    # Get information from input parameter list
    x0_bg,M_4d,H_4d,R_4d = params
    x0_bg = np.matrix(x0_bg).flatten().T

    if verbose:
      print('M_4d = ')
      print(M_4d)
      print('H_4d = ')
      print(H_4d)
      print('R_4d = ')
      print(R_4d)

    # Get model dimension
    xdim = self.xdim

    # Set identity matrix
    I = np.matrix(np.identity(xdim))

    # Get or construct 4D versions of key matrices
    B = np.matrix(self.B)

    # Convert inputs to make sure they are processed as matrices
    Xb = np.matrix(np.atleast_2d(xb_4d))
    Yo = np.matrix(np.atleast_2d(yo_4d))

    if verbose:
      print('Xb = ')
      print(Xb)
      print('Yo = ')
      print(Yo)

    # Track the background estimate from the start of the first loop
    x0_bg = np.matrix(x0_bg)
    if verbose:
      print('x0_bg = ', x0_bg)

    # Get dimensions of input matrices
    Xb_shape = Xb.shape
    Yo_shape = Yo.shape
    M4d_shape = len(M_4d)

    tdim = int(Xb_shape[1])
    ydim = int(Yo_shape[0])

    if verbose:
      print('tdim = ', tdim)
      print('ydim = ', ydim)
      print('Yo = ')
      print(Yo)

    # Set up H_transpose and model equivalent for all times
    Ht = []
    Yb = np.matrix(np.zeros((ydim,tdim)))
    for j in range(tdim):

      # Get linearized h operator for this time
      if verbose:
        print ('j = ', j)

      Ht.append(np.transpose(H_4d[j]))

      # Get model equivalent at observation locations
      if verbose:
        print ('H_4d[j] = ')
        print (H_4d[j])
        print ('Xb[:,j] = ')
        print (Xb[:,j])
     
      Yb[:,j] = np.dot(H_4d[j],Xb[:,j])

      if verbose:
        print('Yb[:,j] = ')
        print(Yb[:,j])

    if verbose:
      print ('Yo = ')
      print (Yo)
      print ('Yb = ')
      print (Yb)

    # Form innovations
    D = np.matrix(np.atleast_2d(Yo - Yb))

    if verbose:
      print ('departures D = ')
      print (D)
      print ('D[0,0] = ')
      print (D[0,0])
#   exit('dep')

    # Pre-compute the Tangent Linear Model M between each time to time t0
    Mt = []
    Rinv = []
    for i in range(tdim):
      Mt.append( deepcopy(np.transpose(M_4d[i])) )
      if verbose:
        print('M[%d] = ' % (i))
        print(M_4d[i])
        print('Mt[%d] = ' % (i))
        print(Mt[i])
        print('TEST: Mt[i]*M[i] = ')
        print(Mt[i]*M_4d[i])
        print('TEST: M[i]*Mt[i] = ')
        print(M_4d[i]*Mt[i])

      # Prepare Rinv list for next section:
      Rinv.append( deepcopy( np.linalg.inv( R_4d[i]) ) )


#   #STEVE: override M to test
#   for i in range(len(M)):
#     M[i] = I
#     Mt[i] = I

#   print('Exiting on purpose...')
#   exit()

    # assign background state on which to optimize (at initial time)
    xb0 = np.matrix(xb_4d[:,0])

    if verbose:
      print ('xb0 = ')
      print (xb0)
    
    # Compute difference with initial background estimate
    # (required for outer loop optimization)
    db0 = x0_bg.A1 - xb0.A1
    db0.reshape((xdim,1))

    if verbose:
      print ('db0 = ')
      print (db0)
#   exit('db0')

    # Compute Jo terms:
    SumMtHtRinvD = np.matrix(np.zeros((xdim,1)))
    for i in reversed(range(tdim)):
    
      if verbose:
        print ('i = ', i)
        print ('D[:,%d] = ' % (i))
        print (D[:,i])

      RinvD = Rinv[i]*D[:,i]

      if verbose:
        print ('RinvD = ')
        print (RinvD)

      HtRinvD = Ht[i]*RinvD

      if verbose:
        print ('HtRinvD = ')
        print (HtRinvD)

      MtHtRinvD = Mt[i]*HtRinvD

      if verbose:
        print ('MtHtRinvD = ')
        print (MtHtRinvD)

      SumMtHtRinvD = SumMtHtRinvD  + MtHtRinvD 

      if verbose:
        print ('SumMtHtRinvD = ')
        print (SumMtHtRinvD)

    if verbose:
      print ('B = ')
      print (B)
      print ('SumMtHtRinvD = ')
      print (SumMtHtRinvD)

    BsumMtHtRinvD = np.matrix(np.dot(B,SumMtHtRinvD))

    if verbose:
      print ('B*SumMtHtRinvD = ')
      print (BsumMtHtRinvD)

    b1 = np.matrix(np.zeros((xdim,1)))

    if verbose:
      print ('db0 = ', db0)
      print ('b1 (init) = ', b1)

    for i in range(xdim):
      b1[i] = db0[i] + BsumMtHtRinvD[i]

    if verbose:
      print ('b1 = ')
      print (b1)

#   exit('Jo exit')

    # Compute Jb terms:
    j = 0
    SumMtHtRinvHM = np.zeros_like(B)
    for i in reversed(range(tdim)):
      if verbose:
        print (' ------------------- i = ', i)
        print ('Ht = ')
        print (Ht[i])
        print ('Mt = ')
        print (Mt[i])

      # Apply obs operator to M
      HM = np.dot(H_4d[i],M_4d[i])
      MtHt = np.transpose(HM)

      # Scale by the obs error
      MtHtRinv = np.dot(MtHt,Rinv[i])

      if verbose:
        print ('MtHt = ')
        print (MtHt)
        print ('MtHtRinv = ')
        print (MtHtRinv)

      # Get symmetric form
      MtHtRinvHM = np.dot(MtHtRinv,HM)

      if verbose:
        print ('MtHtRinvHM = ')
        print (MtHtRinvHM)

      # Sum over time steps
      SumMtHtRinvHM = SumMtHtRinvHM + MtHtRinvHM 

      if verbose:
        print ('SumMtHtRinvHM = ')
        print (SumMtHtRinvHM)

    # Calculate non-symmetric form of 'A' matrix
    B0SumMtHtRinvHM = B*SumMtHtRinvHM
    A = I + B0SumMtHtRinvHM

    #STEVE: consider alternative symmetric formulation to permit use of pcg...

    if verbose:
      print ('4DVar Pre-minimization:')
      print ('A = ')
      print (A)
      print ('b1 = ')
      print (b1)
      print ('xb0 = ')
      print (xb0)

    # Solve for xa (x0):
#   dx0,ierr = sp.sparse.linalg.cg(A, b1, x0=np.zeros_like(b1), tol=1e-05, maxiter=1000)
    # Steve: cg does not converge for large condition numbers of A. bicgstab has been more successful:
    dx0,ierr = sp.sparse.linalg.bicgstab(A,b1,x0=np.zeros_like(b1), tol=1e-05, maxiter=1000)

    dx0 = np.matrix(dx0)
    
    xa = xb0.flatten().T + dx0.flatten().T

    if verbose:
      print ('dx0 = ')
      print (dx0)
      print ('ierr = ')
      print (ierr)
      print('xb0 = ')
      print(xb0)
      print('dx0 = ')
      print(dx0)
      print('xb0.flatten().T = ')
      print(xb0.flatten().T)
      print('dx0.flatten().T = ')
      print(dx0.flatten().T)
      print('xa = ')
      print(xa)

    cn = np.linalg.cond(A)

    if verbose:
      print ('condition number of A = ')
      print (cn)

#   exit('data_assimilation::4dvar::post-cg:: exiting on purpose...')
  
    if (ierr>0):
      print ('Warning in 4DVar conjugate gradient solver, desired tolerance not reached.')
      print ('ierr = ', ierr)
    elif (ierr<0):
      print ('Error in 4DVar conjugate gradient solver. EXITING...')
      print ('ierr = ', ierr)
      exit()

    # Compute KH (based on EKF):
    #K = B_(i+1,i)*H_(i+1)'*( H_(i+1)*B_(i+1,i)*H'_(i+1) + R_(i+1) )^(-1)
#   i=0
#   HM = np.dot(H_4d[i],M_4d[i])
#   MtHt = np.transpose(HM)
#   BMtHt = np.dot(B,MtHt)
#   HMBMtHt = np.dot(HM,BMtHt)
#   S = HMBMtHt + R_4d[i]
#   MBMtHt = np.dot(M_4d[i],BMtHt)
#   Sinv = np.linalg.inv(S)
#   K = np.dot(MBMtHt,Sinv)
#   KH = np.dot(K,H_4d[i])
    K = []
    KH = []

    if verbose:
      print('_4DVarInnerLoop :: -----------------------------------------------------------')
      print('B = ')
      print(B)
      print('R = ')
      print(R_4d[0])
      print('Rinv = ')
      print(Rinv[0])
      print('H = ')
      print(H_4d[0])
      print('Ht = ')
      print(Ht[0])
      print('K = ')
      print(K)
      print('KH = ')
      print(KH)
      print('xb0 = ', xb0)
      print('xa0 = ', xa)
      print('ierr = ', ierr)

#   exit('EXITING on purpose...')

    return np.squeeze(xa).A1,KH


  #-------------------------------------------------------------------------------------------------
  def _4DEnVar(self,Xb_4d,yo_4d,params):
  #-------------------------------------------------------------------------------------------------
  # Use ensemble of states over a time window to estimate temporal correlations
#
# (work in progress)

    verbose = True

    # Get ensemble dimension
    edim = self.edim

    # Format inputs by stacking list of ensemble data as a single state vector for each member
    if verbose:
      print('Xb_4d = ')
      print(Xb_4d)
    Xb = self.convertList2Stacked(Xb_4d)
    if verbose:
      print('Xb = ')
      print(Xb)
    xb = np.mean(Xb,axis=1).T
    if verbose:
      print('xb = ')
      print(xb)

    # Compute forecast perturbations
    Xp = Xb - np.matlib.repmat(xb,1,edim)

    # Get information from input parameter list
    H_4d,R_4d = params

    if verbose:
      print('H_4d = ')
      print(H_4d)
      print('R_4d = ')
      print(R_4d)

    # Update the B, H, and R matrices to correspond to the previous change
    B2 = np.dot(Xp,np.transepose(Xp))
    H2 = self.expandStackedMatrix(H_4d)
    R2 = self.expandStackedMatrix(R_4d)
    das = deepcopy(self)  # Don't modify the original H and R matrices
    das.setB(B2)
    das.setH(H2)
    das.setR(R2)     

    # Call 3DVar on the 4D-augmented system
    xa, KH = das._3DVar(xb.A1,yo_4d)
    Xa = Xp + np.matlib.repmat(xa,1,edim)

    # Select only the appropriate (i.e. final) timestep to output
    return Xa[-xdim:,:], KH[-xdim:,:]


  #-------------------------------------------------------------------------------------------------
  def _4DETKF(self,Xb_4d,yo_4d,params):
  #-------------------------------------------------------------------------------------------------
  # Use ensemble of states over a time window to estimate temporal correlations
#
# (work in progress)
    verbose = True

    # Format inputs by stacking list of ensemble data as a single state vector for each member
    xdim,edim = np.shape(Xb_4d[0])
    Xb = self.convertList2Stacked(Xb_4d)

    # Get information from input parameter list
    H_4d,R_4d = params

    if verbose:
      print('H_4d = ')
      print(H_4d)
      print('R_4d = ')
      print(R_4d)

    # Update the H and R matrices to correspond to the previous change
    H2 = self.expandStackedMatrix(H_4d)
    R2 = self.expandStackedMatrix(R_4d)
    das = deepcopy(self)  # Don't modify the original H and R matrices
    das.setH(H2)
    das.setR(R2)     

    # Call the ETKF on the 4D-augmented system
    Xa, KH = das.ETKF(Xb,yo_4d)

    # Select only the appropriate timestep to output
    return Xa[-xdim:,:], KH[-xdim:,:]

  #-------------------------------------------------------------------------------------------------
  def PF(self,Xb,yo):
  #-------------------------------------------------------------------------------------------------
  # Use an ensemble of states to estimate a multidimensional probability distribution

    verbose = False

    # Make sure inputs are matrix and column vector types
    Xb = np.matrix(Xb)
    yo = np.matrix(yo).flatten().T

    nr,nc = np.shape(Xb)
    xdim = nr
    edim = nc
    ydim = len(yo)
    Hl = self.H
    R = self.R
    Rinv = np.linalg.inv(R)

    # Mintail as machine epsilon
    mintail = np.finfo(float).eps
  
    # Convert background to observation space
    Yb = np.matrix(np.zeros([ydim,edim]))
    for k in range(edim):
      Yb[:,k] = Hl*Xb[:,k]

    # Compute the weights
    # Use a Gaussin for the observation likelihood
    likelihood = np.zeros(edim)
    for k in range(edim):
      likelihood[k] = np.exp(-0.5* (yo-Yb[:,k]).T * Rinv * (yo-Yb[:,k]) )
    likelihood = likelihood + mintail

    # Normalize the weights
    weights = likelihood/np.sum(likelihood)

    # form cumulative distribution
    weight = np.cumsum(weights)

    # Calculate effective sample size
    Neff = 1/np.sum(np.square(weights))

    # Resample the particles with replacement
    addit=1.0/edim
    stt=addit*np.random.rand(1)
    selection_points=np.arange( stt , stt + edim*addit , addit) # Check to make sure ported correctly from matlab
    #DANGER: arange is inconsistent, use linspace instead

    if verbose:
      print('addit = ', addit)
      print('stt = ', stt)
      print('selection_points = ')
      print(selection_points)

    # Set up comb for random selection
    Xa = np.matrix(np.zeros_like(Xb))
    resampling_index = np.zeros(edim)
    j=1; 
    for i in range(edim):
      while selection_points[i] >= weight[j]:
        j=j+1
      resampling_index[i]=j
      Xa[:,i] = Xb[:,j]

    # Apply inflation
#   if (Neff < edim/2):
    if (Neff < edim):
      print('Neff = ', Neff)
      # Apply additive inflation (remove sample mean)
      const=1.0
      stdv_min = 0.01
      # -------------------------------------
      # Correction thanks to Takuma Yoshida:
      # Compute analysis ensemble spread
      stdv = np.std(Xa,axis=1).A
      print('stdv = ', stdv)
      stdv = np.maximum(stdv,stdv_min)
      print('stdv = ', stdv)
      stdv = np.repeat(stdv * const, edim, axis=1)
      print('stdv * const = ', stdv)
      rmat = np.asmatrix(np.random.randn(xdim,edim) * stdv )
      print('rmat = ', rmat)
      print('Xa = ', Xa)
      Xa = Xa + rmat - np.repeat(np.mean(rmat,axis=1),edim,axis=1)
      print('Xa = ', Xa)
      # -------------------------------------
#     exit()
#     rmat=np.random.randn(xdim,edim) * np.matlib.repmat(np.std(Xa,axis=1),1,edim) * const
#     Xa = Xa + rmat - np.matlib.repmat(np.mean(rmat,axis=1),1,edim)

    KH = [0] # dummy output

    return Xa, KH

  #-------------------------------------------------------------------------------------------------
  def HybridGain(self,Xb,yo):
  #-------------------------------------------------------------------------------------------------
  # Use a hybrid method to compute the analysis

    verbose = False

    # Make sure inputs are matrix and column vector types
    Xb = np.matrix(Xb)
    yo = np.matrix(yo).flatten().T

    # Get parameters
    alpha = self.alpha
    nr,nc = np.shape(Xb)
    xdim = nr
    edim = nc
    ydim = len(yo)

    # Compute ETKF analysis
    Xa_ETKF, KH_ETKF = self.ETKF(Xb,yo)
    if verbose:
      print('Xa_ETKF = ')
      print(Xa_ETKF)

    # Compute 3DVar analysis
    xa_ETKF = np.mean(Xa_ETKF,axis=1)

    if verbose:
      print('xa_ETKF = ')
      print(xa_ETKF)

    xa_3DVar, KH_3DVar = self._3DVar(xa_ETKF,yo)

    if verbose:
      print('xa_3DVar = ')
      print(xa_3DVar)

    xa_ETKF = np.matrix(xa_ETKF).flatten().T
    xa_3DVar = np.matrix(xa_3DVar).flatten().T

    # Recover ensemble perturbations
    Xa_ETKF = Xa_ETKF - np.matlib.repmat(xa_ETKF,1,edim)

    # Form hybrid combination to update the mean state
    xa_hybrid = (1-alpha)*xa_ETKF + alpha*xa_3DVar
    Xa = Xa_ETKF + np.matlib.repmat(xa_hybrid,1,edim)

    if verbose:
      print('xa_hybrid = ')
      print(xa_hybrid)

    KH = (1-alpha)*KH_ETKF + alpha*KH_3DVar

    return Xa, KH 


  #-------------------------------------------------------------------------------------------------
  def save(self,outfile):
  #-------------------------------------------------------------------------------------------------
    with open(outfile,'wb') as output:
      pickle.dump(self,output,pickle.HIGHEST_PROTOCOL)

  #-------------------------------------------------------------------------------------------------
  def load(self,infile):
  #-------------------------------------------------------------------------------------------------
    with open(infile,'rb') as input:
      das = pickle.load(input)
    return das
