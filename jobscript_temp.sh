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

echo "new 2*2"
echo "new 2*2" 1>&2
mpirun -map-by ppr:4:socket --bind-to core -n 4 bin/main_new
echo
COMMENTOUT

echo "new 3*3"
echo "new 3*3" 1>&2
mpirun -map-by ppr:9:socket --bind-to core -n 9  bin/main_new
echo

echo "new 6*6"
echo "new 6*6" 1>&2
mpirun -map-by ppr:13:socket --bind-to core -n 36  bin/main_new
echo
