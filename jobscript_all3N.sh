#!/bin/zsh
#SBATCH -p bdw2-mixed
#SBATCH -N 3
#SBATCH -o report/out_%J
#SBATCH -e report/err_%J

echo "new 2*2"
./echo_avetime.sh "mpirun -map-by node --bind-to socket -n 4 -report-bindings bin/main_new"

#echo "old 2*2"
#./echo_avetime.sh "mpirun -map-by ppr:1:node --bind-to core -n 4 -report-bindings bin/main_old"



echo "new 4*4"
./echo_avetime.sh "mpirun -map-by ppr:3:socket --bind-to core -n 16 -report-bindings bin/main_new"

#echo "old 4*4"
#./echo_avetime.sh "mpirun -map-by ppr:6:socket --bind-to core -n 16 -report-bindings bin/main_old"



echo "new 8*8"
./echo_avetime.sh "mpirun -map-by ppr:11:socket --bind-to core -n 64 -report-bindings bin/main_new"

#echo "old 8*8"
#./echo_avetime.sh "mpirun -map-by ppr:22:socket --bind-to core -n 64 -report-bindings bin/main_old"