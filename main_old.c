/*
 * Laplace equation with explicit method
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
/* square region */
#define XSIZE 256
#define YSIZE 256
#define PI 3.1415927
#define NITER 10000
double u[XSIZE + 2][YSIZE + 2], uu[XSIZE + 2][YSIZE + 2];

double **u_old;
double **u_new;

double time1, time2;
void lap_solve(MPI_Comm);
int myid, numprocs;
int namelen;
char processor_name[MPI_MAX_PROCESSOR_NAME];
int xsize;

void initialize()
{
	int x, y;
	/* 初期値を設定*/
	for (x = 1; x < XSIZE + 1; x++)
		for (y = 1; y < YSIZE + 1; y++)
			u[x][y] = sin((x - 1.0) / XSIZE * PI) +
					  cos((y - 1.0) / YSIZE * PI);
	/* 境界をゼロクリア*/
	for (x = 0; x < XSIZE + 2; x++)
	{
		u[x][0] = u[x][YSIZE + 1] = 0.0;
		uu[x][0] = uu[x][YSIZE + 1] = 0.0;
	}
	for (y = 0; y < YSIZE + 2; y++)
	{
		u[0][y] = u[XSIZE + 1][y] = 0.0;
		uu[0][y] = uu[XSIZE + 1][y] = 0.0;
	}
}

#define TAG_1 100
#define TAG_2 101
#ifndef FALSE
#define FALSE 0
#endif
void lap_solve(MPI_Comm comm)
{
	int x, y, k;
	double sum;
	double t_sum;
	int x_start, x_end;
	MPI_Request req1, req2;
	MPI_Status status1, status2;
	MPI_Comm comm1d;
	int down, up;
	int periods[1] = {FALSE};
	/*
	 * Create one dimensional cartesian topology with
	 * nonperiodical boundary
	 */
	MPI_Cart_create(comm, 1, &numprocs, periods, FALSE, &comm1d);
	/* calculate process ranks for 'down' and 'up' */
	MPI_Cart_shift(comm1d, 0, 1, &down, &up);
	x_start = 1 + xsize * myid;
	x_end = 1 + xsize * (myid + 1);
	for (k = 0; k < NITER; k++)
	{
		/* old <- new */
		for (x = x_start; x < x_end; x++)
			for (y = 1; y < YSIZE + 1; y++)
				uu[x][y] = u[x][y];
		/* recv from down */
		MPI_Irecv(&uu[x_start - 1][1], YSIZE, MPI_DOUBLE,
				  down, TAG_1, comm1d, &req1);
		/* recv from up */
		MPI_Irecv(&uu[x_end][1], YSIZE, MPI_DOUBLE,
				  up, TAG_2, comm1d, &req2);
		/* send to down */
		MPI_Send(&u[x_start][1], YSIZE, MPI_DOUBLE,
				 down, TAG_2, comm1d);
		/* send to up */
		MPI_Send(&u[x_end - 1][1], YSIZE, MPI_DOUBLE,
				 up, TAG_1, comm1d);
		MPI_Wait(&req1, &status1);
		MPI_Wait(&req2, &status2);
		/* update */
		for (x = x_start; x < x_end; x++)
			for (y = 1; y < YSIZE + 1; y++)
				u[x][y] = .25 * (uu[x - 1][y] + uu[x + 1][y] +
								 uu[x][y - 1] + uu[x][y + 1]);
	}
	/* check sum */
	sum = 0.0;
	for (x = x_start; x < x_end; x++)
		for (y = 1; y < YSIZE + 1; y++)
			sum += uu[x][y] - u[x][y];
	MPI_Reduce(&sum, &t_sum, 1, MPI_DOUBLE, MPI_SUM, 0, comm1d);
	if (myid == 0)
		fprintf(stderr,"sum = %g\n", t_sum);
	MPI_Comm_free(&comm1d);
}

int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Get_processor_name(processor_name, &namelen);
	fprintf(stderr, "Process %d on %s\n", myid, processor_name);
	xsize = XSIZE / numprocs;
	if ((XSIZE % numprocs) != 0)
		MPI_Abort(MPI_COMM_WORLD, 1);
	initialize();
	MPI_Barrier(MPI_COMM_WORLD);
	time1 = MPI_Wtime();
	lap_solve(MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	time2 = MPI_Wtime();
	if (myid == 0)
	{
		fprintf(stderr,"time = %g\n", time2 - time1);
		printf("%g", time2 - time1);
	}
		
	MPI_Finalize();
	return (0);
}