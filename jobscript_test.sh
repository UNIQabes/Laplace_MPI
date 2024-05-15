#!/bin/zsh
#SBATCH -p bdw2-mixed
#SBATCH -N 3
#SBATCH -o report/out_%J
#SBATCH -e report/err_%J

echo "new 2*2"
mpirun -map-by ppr:4:socket --bind-to core -n 4 -report-bindings bin/main_new

echo "old 2*2"
mpirun -map-by ppr:4:socket --bind-to core -n 4 -report-bindings bin/main_old



echo "new 3*3"
mpirun -map-by ppr:9:socket --bind-to core -n 9 -report-bindings bin/main_new

echo "old 3*3"
mpirun -map-by ppr:9:socket --bind-to core -n 9 -report-bindings bin/main_old




echo "new 4*4"
mpirun -map-by ppr:8:socket --bind-to core -n 16 -report-bindings bin/main_new

echo "old 4*4"
mpirun -map-by ppr:8:socket --bind-to core -n 16 -report-bindings bin/main_old




echo "new 5*5"
mpirun -map-by ppr:13:socket --bind-to core -n 25 -report-bindings bin/main_new

echo "old 5*5"
mpirun -map-by ppr:13:socket --bind-to core -n 25 -report-bindings bin/main_old




echo "new 6*6"
mpirun -map-by ppr:13:socket --bind-to core -n 36 -report-bindings bin/main_new

echo "old 6*6"
mpirun -map-by ppr:13:socket --bind-to core -n 36 -report-bindings bin/main_old




echo "new 7*7"
mpirun -map-by ppr:13:socket --bind-to core -n 49 -report-bindings bin/main_new

echo "old 7*7"
mpirun -map-by ppr:13:socket --bind-to core -n 49 -report-bindings bin/main_old





echo "new 8*8"
mpirun -map-by ppr:13:socket --bind-to core -n 64 -report-bindings bin/main_new

echo "old 8*8"
mpirun -map-by ppr:13:socket --bind-to core -n 64 -report-bindings bin/main_old