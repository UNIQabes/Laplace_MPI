bin/main_old : main_old.c
	mpicc -I /opt/homebrew/opt/libomp/include -L /opt/homebrew/opt/libomp/lib -Xpreprocessor -fopenmp -lomp main_old.c  -o bin/main_old

bin/main_new : main_new.c
	mpicc -I /opt/homebrew/opt/libomp/include -L /opt/homebrew/opt/libomp/lib -Xpreprocessor -fopenmp -lomp main_new.c  -o bin/main_new