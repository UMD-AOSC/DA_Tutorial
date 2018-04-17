#!/bin/bash

set -e

wdir_base="/lustre/tyoshida/shrt/exec"
modeldir="/lustre/tyoshida/prgm/MAOOAM/fortran"

# preparation
word=m`/lustre/tyoshida/repos/python/pythonpath/oneliner/serial`
wdir="${wdir_base}/${word}"

cd ${modeldir}
git commit -a --allow-empty -m "MAOOAM/fortran exec.sh auto commit: experiment ${word}"

echo "preparing files at ${wdir}"
rm -rf ${wdir}
mkdir -p ${wdir}
cd ${wdir}

cp -f ${modeldir}/maooam .
cp -f ${modeldir}/run_true_and_tlm .
cp -f ${modeldir}/*.nml .

echo "#!/bin/bash"                  > tmp.sh
echo "#SBATCH -n 20"                >> tmp.sh
echo "#SBATCH -t 00:15:00"         >> tmp.sh
echo "#SBATCH -J ${word}"          >> tmp.sh
echo "set -e"                      >> tmp.sh
echo "export OMP_NUM_THREADS=20"    >> tmp.sh
# echo "./maooam"                    >> tmp.sh
echo "time ./run_true_and_tlm"          >> tmp.sh

sbatch tmp.sh
