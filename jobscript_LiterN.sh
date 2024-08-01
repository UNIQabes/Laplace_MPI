#!/bin/zsh
#SBATCH -p bdw2-mixed
#SBATCH -N 2
#SBATCH -o report/%J_out_LiterN.txt
#SBATCH -e report/%J_err_LiterN.txt

: << COMMENTOUT
echo "old 4*4"
echo "old 4*4" 1>&2
mpirun -map-by ppr:8:socket --bind-to core -n 16  bin/main_old
echo
COMMENTOUT

echo "LIter2 2*2"
echo "LIter2 2*2" 1>&2
mpirun -map-by ppr:4:socket --bind-to core -n 4 bin/main_LIter_2
echo

echo "LIter3 2*2"
echo "LIter3 2*2" 1>&2
mpirun -map-by ppr:4:socket --bind-to core -n 4 bin/main_LIter_3
echo

echo "LIter4 2*2"
echo "LIter4 2*2" 1>&2
mpirun -map-by ppr:4:socket --bind-to core -n 4 bin/main_LIter_4
echo

echo "LIter5 2*2"
echo "LIter5 2*2" 1>&2
mpirun -map-by ppr:4:socket --bind-to core -n 4 bin/main_LIter_5
echo


: << COMMENTOUT
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
COMMENTOUT


: << COMMENTOUT
echo "LIter 8*8"
echo "LIter 8*8" 1>&2
mpirun -map-by ppr:13:socket --bind-to core -n 64  bin/main_LIter
echo
COMMENTOUT