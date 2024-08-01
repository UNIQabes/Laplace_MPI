LOCALITER=3

MPICC = mpicc
CFLAGS = -std=gnu99 
INCDIR = -I/home/uniqabes/MyCommand/include
LIBDIR  = -L/home/uniqabes/MyCommand/lib
LIB = -ltlog -lm

bin/main_LIter_%: main_LIter.c
	${MPICC} -DLOCALITER=$* ${CFLAGS} ${INCDIR} ${LIBDIR}  $< -o $@ ${LIB} 


bin/%: %.c
	${MPICC} ${CFLAGS} ${INCDIR} ${LIBDIR} $< -o $@ ${LIB} 


all: bin/main_LIter bin/main_new bin/main_old

NLIter: bin/main_LIter_2 bin/main_LIter_3 bin/main_LIter_4 bin/main_LIter_5

.PHONY:all