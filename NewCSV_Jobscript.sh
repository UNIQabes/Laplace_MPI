#!/bin/zsh
#SBATCH -p bdw2-mixed
#SBATCH -N 3
#SBATCH -o report/%J_out_NewCSV.txt
#SBATCH -e report/%J_err_NewCSV.txt

echo "New 2*2"
echo "New 2*2" 1>&2
./make_csv.sh "mpirun -map-by ppr:4:socket --bind-to core -n 4 bin/main_new" "all,Cal,BufToAry,HaloEx,AryToBuf"
echo

echo "New 3*3"
echo "New 3*3" 1>&2
./make_csv.sh "mpirun -map-by ppr:9:socket --bind-to core -n 9  bin/main_new" "all,Cal,BufToAry,HaloEx,AryToBuf"
echo

echo "New 4*4"
echo "New 4*4" 1>&2
./make_csv.sh "mpirun -map-by ppr:8:socket --bind-to core -n 16  bin/main_new" "all,Cal,BufToAry,HaloEx,AryToBuf"
echo

echo "New 5*5"
echo "New 5*5" 1>&2
./make_csv.sh "mpirun -map-by ppr:13:socket --bind-to core -n 25  bin/main_new" "all,Cal,BufToAry,HaloEx,AryToBuf"
echo

echo "New 6*6"
echo "New 6*6" 1>&2
./make_csv.sh "mpirun -map-by ppr:13:socket --bind-to core -n 36  bin/main_new" "all,Cal,BufToAry,HaloEx,AryToBuf"
echo

echo "New 7*7"
echo "New 7*7" 1>&2
./make_csv.sh "mpirun -map-by ppr:13:socket --bind-to core -n 49  bin/main_new" "all,Cal,BufToAry,HaloEx,AryToBuf"
echo

echo "New 8*8"
echo "New 8*8" 1>&2
./make_csv.sh "mpirun -map-by ppr:13:socket --bind-to core -n 64  bin/main_new" "all,Cal,BufToAry,HaloEx,AryToBuf"
echo