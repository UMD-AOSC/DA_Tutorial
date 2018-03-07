# Modular arbitrary-order ocean-atmosphere model: MAOOAM -- Python implementation #

## About ##

(c) 2017 Maxime Tondeur and Jonathan Demaeyer

See [LICENSE.txt](./LICENSE.txt) for license information.

This software is provided as supplementary material with:

* De Cruz, L., Demaeyer, J. and Vannitsem, S.: The Modular Arbitrary-Order
Ocean-Atmosphere Model: MAOOAM v1.0, Geosci. Model Dev., 9, 2793-2808,
[doi:10.5194/gmd-9-2793-2016](http://dx.doi.org/10.5194/gmd-9-2793-2016), 2016.

**Please cite this article if you use (a part of) this software for a
publication.**

The authors would appreciate it if you could also send a reprint of
your paper to <lesley.decruz@meteo.be>, <jonathan.demaeyer@meteo.be> and
<svn@meteo.be>. 

Consult the MAOOAM [code repository](http://www.github.com/Climdyn/MAOOAM)
for updates, and [our website](http://climdyn.meteo.be) for additional
resources.

This python version of the code is available in the [Data Assimilation module DAPPER](https://github.com/nansencenter/DAPPER) .

------------------------------------------------------------------------

## Installation ##

The program can be run with python 2.7 or 3.5. Please note that python 3.5 is needed by the Data Assimilation module DAPPER.
Optionally, F2py is also needed to compile the optimized fortran part.

To install, unpack the archive in a folder or clone with git:

```
git clone https://github.com/Climdyn/MAOOAM.git
```

and run f2py (optional):

```
f2py -c sparse_mult.pyf sparse_mult.f90
```

------------------------------------------------------------------------

##  Description of the files ##

The model tendencies are represented through a tensor called aotensor which
includes all the coefficients. This tensor is computed once at the program
initialization.

* maooam.py : Main program.
* params_maooam.py : A module for the parameters of the model (dimensional, integration and physical).
* ic.py : Initial conditions file.
* ic_def.py : A module that generate the initial conditions if it does not exist.
* inprod_analytic.py : A module that compute the inner product needed for the tensor computation.
* aotensor.py : A module that compute the tensor.
* integrator.py : A module that compute one step of the model and RK2 integration.
* sparse_mult.f90 and sparse_mult.pyf : fortran and f2py files to call fortran module in python.

A documentation is available [here](./doc/build/html/index.html) (html) and [here](./doc/build/latex/MAOOAM.pdf) (pdf).
 
------------------------------------------------------------------------

## Usage ##

The user first has to fill the params_maooam.py according to their needs. See the documentation for more information about this file.
Some examples related to already published articles are available in the params folder.

Finally, the ic.py file specifying the initial condition should be defined. To
obtain an example of this configuration file corresponding to the model you
have previously defined, simply delete the current ic.py file (if it exists)
and run the program:

```
python maooam.py
```

It will generate a new one and start with the 0 initial condition. If you want another 
initial condition, stop the program, fill the newly generated file and restart:

```
python maooam.py
```

The code will generate a file evol_field.dat containing the recorded time evolution of the variables.

------------------------------------------------------------------------

## Implementation notes ##

As the system of differential equations is at most bilinear in y[j] (j=1..n), y
being the array of variables, it can be expressed as a tensor contraction
(written using Einstein convention, i.e. indices that occur twice on one side
of an equation are summed over):

    dy  / dt =  T        y   y      (y  == 1)
      i          i,j,k    j   k       0

The tensor T that encodes the differential equations is composed so that:

* T[i][j][k] contains the contribution of dy[i]/dt proportional to y[j]*y[k].
* Furthermore, y[0] is always equal to 1, so that T[i][0][0] is the constant
contribution to var dy[i]/dt.
* T[i][j][0] + T[i][0][j] is the contribution to  dy[i]/dt which is linear in
y[j].

Ideally, the tensor is composed as an upper triangular matrix 
(in the last two coordinates).

The tensor for this model is composed in the aotensor_def module and uses the
inner products defined in the inprod_analytic module.


------------------------------------------------------------------------

## Final Remarks ##

The code has been translated in Python to be integrated within the [Data Assimilation module DAPPER](https://github.com/nansencenter/DAPPER) .

