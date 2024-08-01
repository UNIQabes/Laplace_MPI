#!/bin/zsh
#SBATCH -p bdw2-mixed
#SBATCH -N 2
#SBATCH -o report/%J_out_LiterN.txt
#SBATCH -e report/%J_err_LiterN.txt



echo "LIter2 2*2"
echo "LIter2 2*2" 1>&2
./make_csv.sh "mpirun -map-by ppr:4:socket --bind-to core -n 4 bin/main_LIter_2" "All,Cal,BufToAry,HaloEx,AryToBuf,Rest"
echo

echo "LIter3 2*2"
echo "LIter3 2*2" 1>&2
./make_csv.sh "mpirun -map-by ppr:4:socket --bind-to core -n 4 bin/main_LIter_3" "All,Cal,BufToAry,HaloEx,AryToBuf,Rest"
echo

echo "LIter4 2*2"
echo "LIter4 2*2" 1>&2
./make_csv.sh "mpirun -map-by ppr:4:socket --bind-to core -n 4 bin/main_LIter_4" "All,Cal,BufToAry,HaloEx,AryToBuf,Rest"
echo

echo "LIter5 2*2"
echo "LIter5 2*2" 1>&2
./make_csv.sh "mpirun -map-by ppr:4:socket --bind-to core -n 4 bin/main_LIter_5" "All,Cal,BufToAry,HaloEx,AryToBuf,Rest"
echo




echo "LIter2 7*7"
echo "LIter2 7*7" 1>&2
./make_csv.sh "mpirun -map-by ppr:13:socket --bind-to core -n 49  bin/main_LIter_2" "All,Cal,BufToAry,HaloEx,AryToBuf,Rest"
echo

echo "LIter3 7*7"
echo "LIter3 7*7" 1>&2
./make_csv.sh "mpirun -map-by ppr:13:socket --bind-to core -n 49  bin/main_LIter_3" "All,Cal,BufToAry,HaloEx,AryToBuf,Rest"
echo

echo "LIter4 7*7"
echo "LIter4 7*7" 1>&2
./make_csv.sh "mpirun -map-by ppr:13:socket --bind-to core -n 49  bin/main_LIter_4" "All,Cal,BufToAry,HaloEx,AryToBuf,Rest"
echo

echo "LIter5 7*7"
echo "LIter5 7*7" 1>&2
./make_csv.sh "mpirun -map-by ppr:13:socket --bind-to core -n 49  bin/main_LIter_5" "All,Cal,BufToAry,HaloEx,AryToBuf,Rest"
echo


