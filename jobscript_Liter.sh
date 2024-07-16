#!/bin/zsh
#SBATCH -p bdw2-mixed
#SBATCH -N 3
#SBATCH -o report/%J_out_test.txt
#SBATCH -e report/%J_err_test.txt

: << COMMENTOUT
echo "old 4*4"
echo "old 4*4" 1>&2
mpirun -map-by ppr:8:socket --bind-to core -n 16  bin/main_old
echo
COMMENTOUT

echo "NIter 2*2"
echo "NIter 2*2" 1>&2
mpirun -map-by ppr:4:socket --bind-to core -n 4 bin/main_LIter
echo


echo "NIter 3*3"
echo "NIter 3*3" 1>&2
mpirun -map-by ppr:9:socket --bind-to core -n 9  bin/main_LIter
echo

echo "NIter 4*4"
echo "NIter 4*4" 1>&2
mpirun -map-by ppr:8:socket --bind-to core -n 16  bin/main_LIter
echo

echo "NIter 5*5"
echo "NIter 5*5" 1>&2
mpirun -map-by ppr:13:socket --bind-to core -n 25  bin/main_LIter
echo

echo "NIter 6*6"
echo "NIter 6*6" 1>&2
mpirun -map-by ppr:13:socket --bind-to core -n 36  bin/main_LIter
echo

echo "NIter 7*7"
echo "NIter 7*7" 1>&2
mpirun -map-by ppr:13:socket --bind-to core -n 49  bin/main_LIter
echo

echo "NIter 8*8"
echo "NIter 8*8" 1>&2
mpirun -map-by ppr:13:socket --bind-to core -n 64  bin/main_LIter
echo
