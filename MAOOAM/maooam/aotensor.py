"""
    Tensor computation module
    =========================

    The equation tensor for the coupled ocean-atmosphere model
    with temperature which allows for an extensible set of modes
    in the ocean and in the atmosphere.

    .. note :: These are calculated using the analytical expressions from De \
    Cruz, L., Demaeyer, J. and Vannitsem, S.: A modular arbitrary-order \
    ocean-atmosphere model: MAOOAM v1.0, Geosci. Model Dev. Discuss.\
            And the `Fortran Code <https://github.com/Climdyn/MAOOAM>`_

    .. note :: The python code is available here : \
    `aotensor.py <../_modules/aotensor.html>`_ .

    :Example:

    >>> aotensor, Li, Lj, Lk, Lv =aotensor.init_aotensor()

    Help Functions
    -------------------
    There are ndim coordinates that correspond to 4 physical quantities.
    These functions help to have the i-th coordinates of a quantity.

    * psi(i) -> i
    * theta(i) -> i + natm
    * A(i) -> i + 2*natm
    * T(i) -> i + 2*natm + noc
    * kdelta(i,j) -> (i==j)

    Global variables
    -------------------

    * real_eps = 2.2204460492503131e-16
    * t=np.zeros( ((ndim+1),(ndim+1),(ndim+1)),dtype=float)

    Dependencies
    ------------------------

    >>> from params_maooam import *
    >>> from inprod_analytic import *
    >>> from scipy.sparse import dok_matrix
    >>> from scipy.sparse import csr_matrix
    >>> import os

    Functions
    ------------------------

    * compute_aotensor
    * coeff(i, j, k, v)
    * simplify
    * init_aotensor


"""

import numpy as np
from scipy.sparse import dok_matrix
from scipy.sparse import csr_matrix

from params_maooam import *
from .inprod_analytic import init_inprod
from .inprod_analytic import atmos
from .inprod_analytic import ocean

real_eps = 2.2204460492503131e-16

t = np.zeros(((ndim+1), (ndim+1), (ndim+1)), dtype=float)


def psi(i):
    return i


def theta(i):
    return i + natm


def A(i):
    return i + 2*natm


def T(i):
    return i + 2*natm + noc


def kdelta(i, j):

    if i == j:
        return 1

    else:
        return 0


def compute_aotensor():
    """
        Computes the three-dimensional tensor t

        Takes the inner products of inprod_analytic \
        and computes the tensor

        :param t: tensor t is a global variable of aotensor
        :type t: array((37,37,37),float)
        :return: change the global tensor
        :rtype: void

        :Example:

        >>> compute_aotensor()

        .. warning:: Needs the global variable aotensor and the global inner \
        products to be initialized.

    """

    t[theta(1), 0, 0] = (Cpa / (1 - atmos.a[0, 0] * sig0))

    for i in range(1, natm+1):
        for j in range(1, natm+1):

            t[psi(i), psi(j), 0] = -((atmos.c[(i-1), (j-1)] * betp) /
                                     atmos.a[(i-1), (i-1)]) - (kd * kdelta((i-1), (j-1))) / 2

            t[theta(i), psi(j), 0] = ((atmos.a[(i-1), (j-1)] * kd * sig0) /
                                      (-2 + 2*atmos.a[(i-1), (i-1)]*sig0))

            t[psi(i), theta(j), 0] = (kd * kdelta((i-1), (j-1))) / 2

            t[theta(i), theta(j), 0] = (-((sig0 * (2. * atmos.c[(i-1), (j-1)] \
            * betp + atmos.a[(i-1), (j-1)] * (kd + 4. * kdp))))\
            + 2.*(LSBpa + sc*Lpa) * kdelta((i-1), (j-1))) \
            / (-2. + 2.*atmos.a[(i-1), (i-1)]*sig0)

            for k in range(1, natm+1):

                t[psi(i), psi(j), psi(k)] = -(atmos.b[(i-1), (j-1), (k-1)] \
                / atmos.a[(i-1), (i-1)])

                t[psi(i), theta(j), theta(k)] = -((atmos.b[(i-1), (j-1), (k-1)] \
                / atmos.a[(i-1), (i-1)]))

                t[theta(i), psi(j), theta(k)] = (atmos.g[(i-1), (j-1), (k-1)] \
                - atmos.b[(i-1), (j-1), (k-1)]*sig0) / \
                (-1 + atmos.a[(i-1), (i-1)]*sig0)

                t[theta(i), theta(j), psi(k)] = (atmos.b[(i-1), (j-1), (k-1)] \
                * sig0) / (1 - atmos.a[(i-1), (i-1)] * sig0)

        for j in range(1, noc+1):

            t[psi(i), A(j), 0] = kd * atmos.d[(i-1), (j-1)] / \
            (2 * atmos.a[(i-1), (i-1)])

            t[theta(i), A(j), 0] = kd * (atmos.d[(i-1), (j-1)] * sig0) \
            / (2 - 2 * atmos.a[(i-1), (i-1)] * sig0)

            t[theta(i), T(j), 0] = atmos.s[(i-1), (j-1)] * (2 * LSBpo + Lpa) \
            / (2 - 2 * atmos.a[(i-1), (i-1)] * sig0)

    for i in range(1, noc+1):
        for j in range(1, natm+1):

            t[A(i), psi(j), 0] = ocean.K[(i-1), (j-1)] * dp \
            / (ocean.M[(i-1), (i-1)] + G)

            t[A(i), theta(j), 0] = -(ocean.K[(i-1), (j-1)]) * dp \
            / (ocean.M[(i-1), (i-1)] + G)

        for j in range(1, noc+1):

            t[A(i), A(j), 0] = -((ocean.N[(i-1), (j-1)] * betp + \
                ocean.M[(i-1), (i-1)] * (rp + dp) * kdelta((i-1), (j-1)))) \
            / (ocean.M[(i-1), (i-1)] + G)

            for k in range(1, noc+1):

                t[A(i), A(j), A(k)] = -(ocean.C[(i-1), (j-1), (k-1)]) \
                / (ocean.M[(i-1), (i-1)] + G)

    for i in range(1, noc+1):

        t[T(i), 0, 0] = Cpo * ocean.W[(i-1), 0]

        for j in range(1, natm+1):

            t[T(i), theta(j), 0] = ocean.W[(i-1), (j-1)] * (2 * sc * Lpo + sbpa)

        for j in range(1, noc+1):

            t[T(i), T(j), 0] = -((Lpo + sbpo)) * kdelta((i-1), (j-1))

            for k in range(1, noc+1):
                t[T(i), A(j), T(k)] = -(ocean.O[(i-1), (j-1), (k-1)])


def coeff(i, j, k, v):
    """
        Affects v for :math:`t_{i,j,k}` making that tensor[i] upper triangular.
          Used in compute_aotensor.

        :param i: first coordinates
        :type i: int in [1,37]
        :param j: second coodinates
        :type j: int in [1,37]
        :param k: third coordinates
        :type k: int in [1,37]
        :param v: value
        :type v: float
        :return: change the global tensor
        :rtype: void

        :Example:

        >>> coeff(i, j, k, v)
     """

    if(abs(v) >= real_eps):

        if j <= k:
            t[i, j, k] = v

        else:
            t[i, k, j] = v


def simplify():
    """
        Make sure that tensor[i] is upper triangular.
        To do after compute_aotensor().

        :param t: tensor t is a global variable of aotensor
        :type t: array((ndim+1,ndim+1,ndim+1),float)
        :return: change the global tensor
        :rtype: void

        :Example:

        >>> simplify()
     """

    for i in range(1, ndim+1):
        for j in range(0, ndim+1):
            for k in range(0, j):

                if t[i, j, k] != 0:
                    t[i, k, j] = t[i, j, k]+t[i, k, j]
                    t[i, j, k] = 0.


def init_aotensor():
    """
        Initialize the tensor.
 
        :return: aotensor, Li, Lj, Lk, Lv
	:type aotensor: (int i, int j, int k, float v) in list
        :type Li: int in list
        :type Lj: int in list
        :type Lk: int in list
        :type Lv: float in list

        :Example:

        >>> aotensor, Li, Lj, Lk, Lv = init_aotensor()
     """
    init_inprod()
    compute_aotensor()
    simplify()
    global tensor_liste
    global tensor
    tensor_liste = []

    for i in range(0, ndim+1):
        Xbis = csr_matrix(t[i])
        X = Xbis
        tensor_liste.append(X)
    tensor = np.array(tensor_liste)

    global aotensor
    aotensor = []
    for i in range(1, ndim+1):
        X = tensor[i].nonzero()
        for m in range(0, tensor[i].nnz):
            j = X[0][m]
            k = X[1][m]
            aotensor.append((i, j, k, tensor[i][j, k]))
    Li = []
    Lj = []
    Lk = []
    Lv = []

    for i in range(0, ndim+1):
        M = t[i]
        Lj1 = M.nonzero()[0]
        Lk1 = M.nonzero()[1]
        for m in range(0, np.count_nonzero(M)):
            j = Lj1[m]
            k = Lk1[m]
            Li.append(i)
            Lj.append(j)
            Lk.append(k)
            Lv.append(M[j, k])
    return aotensor, Li, Lj, Lk, Lv

