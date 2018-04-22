# DA Tutorial with MAOOAM
## Usage
```bash
cd maooam_fortran && make && cd ../ && cp maooam_fortran/step_maooam.so .
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
    * Included in runall_tutorial_1.sh (this uses pure python and is slow)

* Analysis core programs for ETKF:
    * **generate_analysis_3dEns.py**
    * **class_da_system.py**

## Todo and Wishlist (feel free to edit me)
### Todo
* Limit B to projection of some modes
    * [X] BV-like mode separation of slower modes
    * [ ] CLV-based mode separation
        * CLVs are, in practice, difficult to obtain in real applications because it uses future dynamics. But it worth examined from theoretical perspective.
    * [X] Calculation of static B
* [ ] Use of multiple covariances/gains
* [ ] Test, test, test...

### Wishlist
* More plots
    * [X] RMS error and spread (ordinate) vs time (abscissa), separately for each component {atm-psi, atm-theta, ocn-psi, ocn-theta}
    * [ ] CLVs
    * [X] B (raw matrix and its eigenvectors)
* [x] Parallelization of ensemble integration by multiprocessing.pool (**generate_analysis_3dEns.py** line 100)
    * Tested with branch "parallel". Not substantial speedup.
* [X] Speedup by using fortran integration
    * About 100x faster (7a7b78d)
    * Note that {int_params.nml, modeselection.nml, params.nml} are needed. Take care those doesn't diverge from parameters for python-MAOOAM.
* [X] Rough estimation of necessary experiment length
    * Ocean streamfunction has timescale of 1E+5 time units (~ 30 years). Experiments with 1E+6 time units are enough.
    * With time step of 0.1 time units, 1E+6 time units (1E+7 steps) single integration is about 80 secs after speedup.
* [ ] Ens members should start from initial conditions independent from truth
    * **analysis_init.py** line 72
    * **generate_analysis_3dEns.py** lines 24-26
    * I heard this only changes the first several windows. I will do it later.

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
