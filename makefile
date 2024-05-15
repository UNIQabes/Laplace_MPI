bin/main_old : main_old.c
	mpicc -lm main_old.c  -o bin/main_old

bin/main_new : main_new.c
	mpicc -lm main_new.c  -o bin/main_new