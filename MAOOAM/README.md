# DA Tutorial with MAOOAM
## Usage
```bash
sh runall.sh
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

## Tested environment
* gfortran 5.4.0
* Anaconda3-5.1.0
    * Python 3.6.4
    * numpy 1.14.0
    * matplotlib 2.1.2
    * scipy 1.0.0
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
