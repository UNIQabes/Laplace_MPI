MPICC = mpicc
CFLAGS = -lm -std=gnu99


bin/%: %.c
	${MPICC} ${CFLAGS} $< -o $@

.PHONY:all

all: bin/main_LIter bin/main_new bin/main_old