"""
        Initial conditions file
        =========================

        This file defines the initial conditions of the model. \
        To be deleted if the dimensions are changed.

        :Example:

        >>> from ic import X0

        Global variables (state vectors)
        --------------------------------

        * X0 random (non-null) initial conditions.

        Dependencies
        ------------

        >>> import numpy as np

"""
import numpy as np

X0 = np.zeros(36)

# psi variables
X0[0] = 0.02  # typ=A, NX0=0, Ny= 1.0
X0[1] = 0.0  # typ=K, NX0=1.0, Ny= 1.0
X0[2] = 0.0  # typ=L, NX0=1.0, Ny= 1.0
X0[3] = 0.0  # typ=A, NX0=0, Ny= 2.0
X0[4] = 0.0  # typ=K, NX0=1.0, Ny= 2.0
X0[5] = 0.0  # typ=L, NX0=1.0, Ny= 2.0
X0[6] = 0.0  # typ=K, NX0=2.0, Ny= 1.0
X0[7] = 0.0  # typ=L, NX0=2.0, Ny= 1.0
X0[8] = 0.0  # typ=K, NX0=2.0, Ny= 2.0
X0[9] = 0.0  # typ=L, NX0=2.0, Ny= 2.0

# theta variables
X0[10] = 0.0  # typ=A, NX0=0, Ny= 1.0
X0[11] = 0.0  # typ=K, NX0=1.0, Ny= 1.0
X0[12] = 0.0  # typ=L, NX0=1.0, Ny= 1.0
X0[13] = 0.0  # typ=A, NX0=0, Ny= 2.0
X0[14] = 0.0  # typ=K, NX0=1.0, Ny= 2.0
X0[15] = 0.0  # typ=L, NX0=1.0, Ny= 2.0
X0[16] = 0.0  # typ=K, NX0=2.0, Ny= 1.0
X0[17] = 0.0  # typ=L, NX0=2.0, Ny= 1.0
X0[18] = 0.0  # typ=K, NX0=2.0, Ny= 2.0
X0[19] = 0.0  # typ=L, NX0=2.0, Ny= 2.0

# A variables
X0[20] = 0.0  # NX0=0.5, Ny= 1.0
X0[21] = 0.0  # NX0=0.5, Ny= 2.0
X0[22] = 0.0  # NX0=0.5, Ny= 3.0
X0[23] = 0.0  # NX0=0.5, Ny= 4.0
X0[24] = 0.0  # NX0=1.0, Ny= 1.0
X0[25] = 0.0  # NX0=1.0, Ny= 2.0
X0[26] = 0.0  # NX0=1.0, Ny= 3.0
X0[27] = 0.0  # NX0=1.0, Ny= 4.0

# T variables
X0[28] = 0.0  # NX0=0.5, Ny= 1.0
X0[29] = 0.0  # NX0=0.5, Ny= 2.0
X0[30] = 0.0  # NX0=0.5, Ny= 3.0
X0[31] = 0.0  # NX0=0.5, Ny= 4.0
X0[32] = 0.0  # NX0=1.0, Ny= 1.0
X0[33] = 0.0  # NX0=1.0, Ny= 2.0
X0[34] = 0.0  # NX0=1.0, Ny= 3.0
X0[35] = 0.0  # NX0=1.0, Ny= 4.0

