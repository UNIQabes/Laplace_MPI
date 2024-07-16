#!/bin/zsh
#SBATCH -p bdw2-mixed
#SBATCH -N 3
#SBATCH -o report/%J_out_Liter.txt
#SBATCH -e report/%J_err_Liter.txt

: << COMMENTOUT
echo "old 4*4"
echo "old 4*4" 1>&2
mpirun -map-by ppr:8:socket --bind-to core -n 16  bin/main_old
echo
COMMENTOUT

echo "LIter 2*2"
echo "LIter 2*2" 1>&2
mpirun -map-by ppr:4:socket --bind-to core -n 4 bin/main_LIter
echo

echo "LIter 3*3"
echo "LIter 3*3" 1>&2
mpirun -map-by ppr:9:socket --bind-to core -n 9  bin/main_LIter
echo

echo "LIter 4*4"
echo "LIter 4*4" 1>&2
mpirun -map-by ppr:8:socket --bind-to core -n 16  bin/main_LIter
echo

echo "LIter 5*5"
echo "LIter 5*5" 1>&2
mpirun -map-by ppr:13:socket --bind-to core -n 25  bin/main_LIter
echo

echo "LIter 6*6"
echo "LIter 6*6" 1>&2
mpirun -map-by ppr:13:socket --bind-to core -n 36  bin/main_LIter
echo

echo "LIter 7*7"
echo "LIter 7*7" 1>&2
mpirun -map-by ppr:13:socket --bind-to core -n 49  bin/main_LIter
echo

echo "LIter 8*8"
echo "LIter 8*8" 1>&2
mpirun -map-by ppr:13:socket --bind-to core -n 64  bin/main_LIter
echo
