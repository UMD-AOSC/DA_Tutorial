"""
    Integration module
    ======================

    This module actually contains the Heun algorithm routines.

    .. note :: The python code is available here : \
    `integrator.py <../_modules/integrator.html>`_ .

    :Example:

    >>> from integrator import step
    >>> step(y,t,dt)

    Global variables
    -------------------

    * **aotensor** tensor with the format (int i, int j, int k, float v) in list
    * **Li** first list of index of tensor
    * **Lj** second list of index of tensor
    * **Lk** third list of index of tensor
    * **Lv** list of tensor values

    Dependencies
    -------------------

    >>> import numpy as np
    >>> from params_maooam import ndim,f2py
    >>> import aotensor as aotensor_mod
    >>> if f2py:
    >>>     import sparse_mult as mult
    >>>     sparse_mul3_f2py = mult.sparse_mult.sparse_mul3

    Functions
    -------------------

    * sparse_mul3
    * tendencies
    * step
"""
import numpy as np

from params_maooam import ndim,f2py
from . import aotensor as aotensor_mod
if f2py:
    import sparse_mult as mult
    sparse_mul3_fortran = mult.sparse_mult.sparse_mul3

aotensor, Li, Lj, Lk, Lv = aotensor_mod.init_aotensor()


def sparse_mul3_py(arr):
    r"""
    Calculate for each i the sums on j,k of the product

    .. math:: 
        tensor(i,j,k)* arr(j) * arr(k)

    .. note:: 
        Python-only function
    """

    if np.ndim(arr) is 1:
        arr = arr.reshape((1, len(arr)))

    n = arr.shape[0]
    arr = np.hstack((np.array([[1]*n]).reshape(n, 1), arr))
    res = np.zeros((n, ndim+1))

    for (i,j,k,v) in aotensor :
         res[:,i]=res[:,i]+v*arr[:,j]*arr[:,k]

    if res.shape[0] is 1:
        res = res.squeeze()
        return res[1:]

    return res[:, 1:]

if f2py:
    def sparse_mul3_f2py(arr):
        r"""
        Calculate for each i the sums on j,k of the product

        .. math:: 
            tensor(i,j,k)* arr(j) * arr(k)

        .. note:: 
            Use the fortran module sparse-mul3_f2py
        """

        if np.ndim(arr) is 1:
            arr = arr.reshape((1, len(arr)))

        n = arr.shape[0]
        arr = np.hstack((np.array([[1]*n]).reshape(n, 1), arr))
        res = np.zeros((n, ndim+1))

        res = sparse_mul3_fortran(Li, Lj, Lk, Lv, arr)

        if res.shape[0] is 1:
            res = res.squeeze()
            return res[1:]

        return res[:, 1:]

if f2py:
    sparse_mul3=sparse_mul3_f2py
else:
    sparse_mul3=sparse_mul3_py

def tendencies(y):
    """ Calculate the tendencies thanks to the product of the tensor and \
    the vector y"""
    return sparse_mul3(y)


def step(y, t, dt):
    """ RK2 method integration"""

    n = y.shape[0]

    buf_f0 = np.zeros((n, ndim+1))
    buf_f1 = np.zeros((n, ndim+1))
    buf_y1 = np.zeros((n, ndim+1))

    buf_f0 = tendencies(y)
    buf_y1 = y + dt * buf_f0
    buf_f1 = tendencies(buf_y1)

    Y = y + 0.5 * (buf_f0 + buf_f1) * dt

    return Y

if __name__ == "__main__":
    import ic
    X = ic.X0
    print(X)
    X = sparse_mul3(X)
    print(X)

