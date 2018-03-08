"""
    Parameters module
    ======================

    This module defines the parameters for the model.

    .. note :: The python code is available here : \
    `params_maooam.py <../_modules/params_maooam.html>`_ .

    :Example:

    >>> from params_maooam import ndim,natm,noc
    >>> from params_maooam import oms,nboc,ams,nbatm
    >>> from params_maooam import *

    There are three types of parameters :

    * integration parameters : simulation time (transient and effective), \
    time step, writeout and write step time
    * dimensional parameters : dimensions of the truncation of fourier \
    for the atmosphere and the ocean
    * physical parameters : they are used in the tensor for the integration

    Integration parameters
    ----------------------

    .. warning:: Time is adimensional. If t_real is in seconds,
                then t_model = t_real * f_0 where f_0 is the Coriolis
                parameter at 45 degrees latitude ( 1.032e-4 )

    * **t_trans** : the transient simulation time of the model to be on the attractor. The states vectors are not written on evol_field.dat.
    * **t_run** : the running simulation time of the model. The states vectors are written on evol_field.dat every tw.
    * **dt** : the step time.
    * **writeout** : boolean value to decide if the module produces evol_field.dat.
    * **tw** : the step time to write on evol_field.
    * **f2py** : boolean to activate the f2py optimization.

    Dimensional parameters
    ----------------------

    * **oms** and **ams** : the matrices that gives the possible values of the modes Nx and Ny.
    * **nboc** and **natm** : the numbers of oceanic and atmospheric blocs.
    * **natm** and **noc** : the numbers of functions available.
    * **ndim** : the total dimension.


    :Example:

    >>> oms =get_modes(2,4)# ocean mode selection
    >>> ams =get_modes(2,2)# atmosphere mode selection
    >>> nboc,nbatm = 2*4,2*2      # number of blocks
    >>> (natm,noc,ndim)=init_params(nboc,nbatm)
    >>>
    >>> # Oceanic blocs
    >>> #( x block number accounts for half-integer wavenumber e.g 1\
    => 1/2 , 2 => 1, etc...)
    >>> OMS[0,:] = 1,1
    >>> OMS[1,:] = 1,2
    >>> OMS[2,:] = 1,3
    >>> OMS[3,:] = 1,4
    >>> OMS[4,:] = 2,1
    >>> OMS[5,:] = 2,2
    >>> OMS[6,:] = 2,3
    >>> OMS[7,:] = 2,4
    >>> #Atmospheric blocs
    >>> AMS[0,:] = 1,1
    >>> AMS[1,:] = 1,2
    >>> AMS[2,:] = 2,1
    >>> AMS[3,:] = 2,2

    :Typical dimensional parameters:

    * atmosphere 2x,2y ; ocean 2x,4y ; ndim = 36
    * atmosphere 2x,2y ; ocean 4x,4y ; ndim = 52
    * atmosphere 2x,4y ; ocean 2x,4y ; ndim = 56

    Physical parameters
    -------------------

    Some defaut parameters are presented below.
    Some parameters files related to already published article are available in the params folder.

    Scale parameters
    ^^^^^^^^^^^^^^^^

    * **scale = 5.e6**  : characteristic space scale, L*pi
    * **f0 = 1.032e-4**  : Coriolis parameter at 45 degrees latitude
    * **n = 1.5e0**  : aspect ratio (n = 2Ly/Lx ; Lx = 2*pi*L/n; Ly = pi*L)
    * **rra = 6370.e3**  : earth radius
    * **phi0_npi = 0.25e0**  : latitude exprimed in fraction of pi

    Parameters for the ocean
    ^^^^^^^^^^^^^^^^^^^^^^^^^^

    * **gp = 3.1e-2**  : reduced gravity
    * **r = 1.e-8**  : frictional coefficient at the bottom of the ocean
    * **h = 5.e2**  : depth of the water layer of the ocean
    * **d = 1.e-8**  : the coupling parameter (should be divided by f0 to be adim)

    Parameters for the atmosphere
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    * **k = 0.02**  : atmosphere bottom friction coefficient
    * **kp = 0.04**  : atmosphere internal friction coefficient
    * **sig0 = 0.1e0**  : static stability of the atmosphere

    Temperature-related parameters for the ocean
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    * **Go = 2.e8**  : Specific heat capacity of the ocean (50m layer)
    * **Co = 350**  : Constant short-wave radiation of the ocean
    * **To0 = 285.0**  : Stationary solution for the 0-th order ocean temperature

    Temperature-related parameters for the atmosphere
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    * **Ga = 1.e7**  : Specific heat capacity of the atmosphere
    * **Ca = 100.e0**  ; Constant short-wave radiation of the atmosphere
    * **epsa = 0.76e0**  : Emissivity coefficient for the grey-body atmosphere
    * **Ta0 = 270.0**  : Stationary solution for the 0-th order atmospheric temperature

    Other temperature-related parameters/constants
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    * **sc = 1.**  : Ratio of surface to atmosphere temperature
    * **lambda = 20.00**  : Sensible+turbulent heat exchange between oc and atm
    * **rr = 287.e0**  : Gas constant of dry air
    * **sb = 5.6e-8**  : Stefan-Boltzmann constant

    Key values
    ^^^^^^^^^^
    * **k** is the friction coefficient at the bottom of the atmosphere. Typical values are 0.01 or 0.0145 for chaotic regimes.
    * **kp** is the internal friction between the atmosphere layers. kp=2*k
    * **d** is the friction coefficient between the ocean and the atmosphere. Typical values are 6*10^{-8} s^{-1} or 9*10^{-8} s^{-1}.
    * **lambda** is the heat exchange between the ocean and the atmosphere. Typical values are 10 W m^{-2} K^{-1} or 15.06 W m ^{-2} K^{-1}.

    Dependencies
    -------------------

    >>> import numpy as np

    Fonctions
    ---------

    Here are the functions to generate the parameters.
"""
import numpy as np


# -----------------------------------------------------------
# Integration parameters
# -----------------------------------------------------------

t_trans = 10  # transient period (e.g. 1.e7)
t_run = 10  # length of trajectory on the attractor (e.g. 5.e8)
dt = 1.e-2  # the time step
writeout = True  # write out all variables every tw time units
tw = 10.0  # the time step of writing
f2py = False # activate the f2py optimization

# -----------------------------------------------------------
# Number of blocks
# -----------------------------------------------------------

# generate the mode blocks for either ocean or atmosphere
# up to given wavenumbers
# nxmax and nymax are the maximum


def get_modes(nxmax, nymax):
    """ Computes the matrix oms and ams with nxmax and nymax"""
    res = np.zeros((nxmax * nymax, 2))
    i = 0
    for Nx in range(1, nxmax+1):
        for Ny in range(1, nymax+1):
            res[i] = [Nx, Ny]
            i += 1
    return res


def init_params(nboc, nbatm):
    """ Computes the dimensions of the system"""
    natmres = 0
    for i in range(0, nbatm):
        if (ams[i, 0] == 1):
            natmres = natmres + 3
        else:
            natmres = natmres + 2
    s = np.shape(oms)
    nocres = s[0]

    return (natmres, nocres, 2 * (natmres + nocres))

# select (.,.) the maximum value admitted for Nx and Ny
# don't forget to delete ic.py, it will regenerates

oms = get_modes(2, 4)  # ocean mode selection
ams = get_modes(2, 2)  # atmosphere mode selection
nboc, nbatm = 2 * 4, 2 * 2  # number of blocks
(natm, noc, ndim) = init_params(nboc, nbatm)


# noc,natm=8,1  # number of basis functions
# ndim=36  # number of variables

# -----------------------------------------------------------
# Scale parameters for the ocean and the atmosphere
# -----------------------------------------------------------

scale = 5.e6  # the characteristic space scale, L*pi
f0 = 1.032e-4  # Coriolis parameter at 45 degrees latitude
n = 1.5e0  # aspect ratio (n = 2Ly/Lx ; Lx = 2*pi*L/n; Ly = pi*L)
rra = 6370.e3  # earth radius
phi0_npi = 0.25e0  # latitude exprimed in fraction of pi

# Parameters for the ocean
gp = 3.1e-2  # reduced gravity
r = 1.e-8  # frictional coefficient at the bottom of the ocean
h = 5.e2  # depth of the water layer of the ocean
d = 1.e-8  # the coupling parameter (should be divided by f0 to be adim)

# Parameters for the atmosphere
k = 0.02  # atmosphere bottom friction coefficient
kp = 0.04  # atmosphere internal friction coefficient
sig0 = 0.1e0  # static stability of the atmosphere

# Temperature-related parameters for the ocean
Go = 2.e8  # Specific heat capacity of the ocean (50m layer)
Co = 350  # Constant short-wave radiation of the ocean
To0 = 285.0  # Stationary solution for the 0-th order ocean temperature

# Temperature-related parameters for the atmosphere
Ga = 1.e7  # Specific heat capacity of the atmosphere
Ca = 100.e0  # Constant short-wave radiation of the atmosphere
epsa = 0.76e0  # Emissivity coefficient for the grey-body atmosphere
Ta0 = 270.0  # Stationary solution for the 0-th order atmospheric temperature

# Other temperature-related parameters/constants
sc = 1.  # Ratio of surface to atmosphere temperature
lambdaa = 20.00  # Sensible+turbulent heat exchange between oc and atm
rr = 287.e0  # Gas constant of dry air
sb = 5.6e-8  # Stefan-Boltzmann constant

# -----------------------------------------------------------
# Some general parameters (Domain, beta, gamma, coupling
# -----------------------------------------------------------

pi = np.arccos(-1.e0)
L = scale / pi
phi0 = phi0_npi * pi
LR = np.sqrt(gp * h) / f0
G = -L**2 / LR**2
betp = L / rra * np.cos(phi0) / np.sin(phi0)
rp = r / f0
dp = d / f0
kd = k * 2
kdp = kp

# -----------------------------------------------------------
# Derived Quantities
# -----------------------------------------------------------

Cpo = Co / (Go * f0) * rr / (f0**2 * L**2)

Lpo = lambdaa / (Go * f0)

# Cpa acts on psi1-psi3, not on theta :
Cpa = Ca / (Ga * f0) * rr / (f0**2 * L**2) / 2

Lpa = lambdaa / (Ga * f0)

# long wave radiation lost by ocean to atmosphere space :
sbpo = 4 * sb * To0**3 / (Go * f0)

# long wave radiation from atmosphere absorbed by ocean
sbpa = 8 * epsa * sb * Ta0**3 / (Go * f0)

# long wave radiation from ocean absorbed by atmosphere
LSBpo = 2 * epsa * sb * To0**3 / (Ga * f0)

# long wave radiation lost by atmosphere to space & ocean
LSBpa = 8 * epsa * sb * Ta0**3 / (Ga * f0)

