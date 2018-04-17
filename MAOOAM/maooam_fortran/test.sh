set -e
rm -f *.pdf *.dat
make clean
make
./maooam
time python call_ft.py
