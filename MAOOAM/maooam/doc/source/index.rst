Modular arbitrary-order ocean-atmosphere model: MAOOAM -- Python implementation
=======================================================================================

Synopsis
-------------------

This repository provides the code of the model MAOOAM in python. It is a low-order ocean-atmosphere model with an arbitrary expansion of the Fourier modes of temperatures and streamfunctions.
The code in Python is a translation of the Fortran code available in the main Git repository : https://github.com/Climdyn/MAOOAM.

Motivation
-------------------

The code has been translated in Python to be used with the Data Assimilation module DAPPER : https://github.com/nansencenter/DAPPER .

Installation
-------------------

The program can be run with python 2.7 or 3.5. Please note that python 3.5 is needed by the Data Assimilation module DAPPER.
Optionally, F2py is also needed to compile the optimized fortran part.

To install, unpack the archive in a folder or clone with git:
     >>> git clone https://github.com/Climdyn/MAOOAM.git
and run f2py (optional):
     >>> f2py -c sparse_mult.pyf sparse_mult.f90

Getting started
-------------------

The user first has to fill the `params_maooam.py <../html/_modules/params_maooam.html>`_ according to their needs. See `its documentation <../html/rstfiles/params_maooam.html>`_ for more information about this file.
Some examples related to already published articles are available in the params folder.

Finally, the ic.py file specifying the initial condition should be defined. To
obtain an example of this configuration file corresponding to the model you
have previously defined, simply delete the current ic.py file (if it exists)
and run the program:

     >>> ipython maooam.py

It will generate a new one and start with the 0 initial condition. If you want another 
initial condition, stop the program, fill the newly generated file and restart:

     >>> ipython maooam.py

The code will generate a file evol_field.dat containing the recorded time evolution of the variables.

Description of the files 
------------------------

* `maooam.py <../html/_modules/maooam.html>`_ : main program.

* `params_maooam.py <../html/_modules/params_maooam.html>`_ : a module for the parameters of the model (dimensional, integral and physical).

* ic.py : initial conditions file.

* `ic_def.py <../html/_modules/ic_def.html>`_ a module that generate the initial conditions if it does not exist.

* `inprod_analytic.py <../html/_modules/inprod_analytic.html>`_ : a module that compute the inner product needed for the tensor computation.

* `aotensor.py <../html/_modules/aotensor.html>`_ : a module that compute the tensor.

* `integrator.py <../html/_modules/integrator.html>`_: a module that compute one step of the model and RK2 integration.

* sparse_mult.f90 and sparse_mult.pyf : fortran and f2py files to call fortran module in python.

Contents
--------------

.. toctree::
   :maxdepth: 2
   :numbered:


   rstfiles/params_maooam
   rstfiles/ic_def
   rstfiles/ic
   rstfiles/inprod_analytic
   rstfiles/aotensor
   rstfiles/integrator
   rstfiles/maooam


Contributors
-------------------

Maxime Tondeur, Jonathan Demaeyer

License
-------------------

© 2017 Maxime Tondeur and Jonathan Demaeyer

See `LICENSE.txt <../../../LICENSE.txt>`_  for license information.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

