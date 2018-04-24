set -e

cd maooam_fortran && make && cd ../ && cp maooam_fortran/step_maooam.so .

sh clean.sh
sh runall_tutorial_1.sh

sh runall_tutorial_2.sh 3DVar
sh runall_tutorial_3.sh ETKF
sh runall_tutorial_3.sh hybrid

mkdir -p img/3DVar img/ETKF img/hybrid
python plot_error.py 3DVar
python plot_error.py ETKF
python plot_error.py hybrid

python module_plot.py
