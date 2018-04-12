# DA Tutorial with MAOOAM
## Usage
```bash
# prepare truth, freerun, and obs
sh clean.sh
sh runall_tutorial_1.sh
# execute ETKF
sh runall_tutorial_3.sh ETKF
# plot (makes directory img/)
python module_plot.py
python plot_error.py ETKF
```

## Key points to be edited often
* Change exp length
**params_maooam.py** -> t_run

* Change obs network (H operator)
    * Edit both of:
        * **generate_observations.py** (line 46) for generation of observations from truth
        * **analysis_init.py** (line 108) for R and H used in DA
    * Current observation error is very small. It should be enlarged for more nonlinearity
    * **module_obs_network.py** provides static H operator with gridpoint obs of U and T. It can be changed.

* To calculate CLVs:
    * I haven't tried with this python version. Should be asked to other students.

* Analysis core programs for ETKF:
    * **generate_analysis_3dEns.py**
    * **class_da_system.py**

## Todo and Wishlist (feel free to edit me)
### Todo
* Limit B to projection of some modes
    * BV-like mode separation of slower modes
    * CLV-based mode separation
    * Static and dynamic B

### Wishlist
* More plots
    * L2 norm error and spread (ordinate) vs time (abscissa), separately for each component {atm-psi, atm-theta, ocn-psi, ocn-theta}
    * CLVs and B (raw matrix and its eigenvectors)
* Parallelization of ensemble integration by multiprocessing.pool
* Rough estimation of necessary experiment length
* Ens members should start from initial conditions independent from truth
    * **analysis_init.py** line 72
    * **generate_analysis_3dEns.py** lines 24-26

## Tested environment
* Anaconda3-5.1.0
    * Python 3.6.4
    * numpy 1.14.0
    * matplotlib 2.1.2
    * scipy 1.0.0
* Part of the program may fail to run with < Python 3.5
    * subprocess.run in **module_plot.py** is a new feature. You can comment them out and execute mkdir manually.
* On Windows Subsystem for Linux

## Original README from https://github.com/UMD-AOSC/DA_Tutorial
```bash
python -m site --user-site
mkdir -p <result>

# add to MAOOAM directory:
touch __init__.py

# copy:
cp <PATH>/MAOOAM/python/params_maooam.py .
```
