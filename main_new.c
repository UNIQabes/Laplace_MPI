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

#define TAG_1 100
#define TAG_2 101
#ifndef FALSE
#define FALSE 0
#endif

int main(int argc, char *argv[])
{
	const int DIM = 2;

	int namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	int myrank, numprocs;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Get_processor_name(processor_name, &namelen);
	fprintf(stderr, "Process %d on %s\n", myrank, processor_name);

	// 各パラメーターの決定--------------------------------------
	//  x,y軸方向に幾つ分割されているか
	int gridNum_x;
	int gridNum_y;

	// 1グリッドが担当する範囲の大きさ
	int gridSize_x;
	int gridSize_y;

	gridNum_x = (int)(sqrtf(numprocs) + 1);
	for (gridNum_x = (int)(sqrtf(numprocs) + 1); gridNum_x >= 1; gridNum_x--)
	{
		if (numprocs % gridNum_x == 0)
		{
			gridNum_y = numprocs / gridNum_y;
			break;
		}
	}
	gridSize_x = XSIZE / gridNum_x;
	gridSize_y = YSIZE / gridNum_y;

	MPI_Comm comm2d;
	int periods[DIM] = {FALSE, FALSE};
	MPI_Cart_create(MPI_COMM_WORLD, 1, &numprocs, periods, FALSE, &comm2d);

	int coords[DIM];
	int maxDim;
	MPI_Cart_coords(comm2d, myrank, 2, coords);
	int gridPos_x = coords[0];
	int gridPos_y = coords[1];

	// 端のグリッドは、gridsize_x gridsize_yよりも担当領域が小さくなることがある。その大きさを求める
	int xstart, xend, ystart, yend;
	xstart = gridSize_x * gridPos_x;
	xend = gridSize_x * (gridPos_x + 1);
	xend = (xend >= XSIZE) ? XSIZE - 1 : xend;
	ystart = gridSize_y * gridPos_y;
	yend = gridSize_y * (gridPos_y + 1);
	yend = (yend >= YSIZE) ? YSIZE - 1 : yend;
	int localGridSize_x = xend - xstart;
	int localGridSize_y = yend - ystart;

	// uoldとunewのメモリ領域確保 + u_oldの初期化 ----------------------------
	//  u_oldとu_newは自分の担当領域+その隣接領域の値を保持する。
	//  xの範囲は gridSize_x*gridPos_?-1 <= x <= gridSize_x*(gridPos_x+1) yの範囲も同様
	double **u_old = (double **)malloc(sizeof(double *) * (localGridSize_x + 2));
	double **u_new = (double **)malloc(sizeof(double *) * (localGridSize_x + 2));

	for (int localx = 0; localx <= localGridSize_x + 1; localx++)
	{
		u_old[localx] = (double *)malloc(sizeof(double) * (localGridSize_y + 2));
		u_new[localx] = (double *)malloc(sizeof(double) * (localGridSize_y + 2));

		int worldx = localx - 1 + gridSize_x * (gridPos_x);

		for (int localy = 0; localy <= localGridSize_y + 1; localy++)
		{
			int worldy = localy - 1 + gridSize_y * (gridPos_y);
			if (worldx < 0 || worldy < 0 || XSIZE - 1 < worldx || YSIZE - 1 < worldy)
			{
				u_old[localx][localy] = 0;
			}
			else
			{
				u_old[localx][localy] = sin(worldx / XSIZE * PI) + cos(worldy / YSIZE * PI);
			}
		}
	}

	// laplaceを解く-----------------------------------
	MPI_Barrier(comm2d);
	double start = MPI_Wtime();

	int xdown, xup, ydown, yup;
	/* calculate process ranks for 'down' and 'up' */
	MPI_Cart_shift(comm2d, 0, 1, &xdown, &xup);
	MPI_Cart_shift(comm2d, 1, 1, &ydown, &yup);

	MPI_Request req_xup, req_xdown, req_yup, req_ydown;

	// 担当領域の一番端(edge)
	double *yup_edge = (double *)malloc(sizeof(double) * localGridSize_x);
	double *ydown_edge = (double *)malloc(sizeof(double) * localGridSize_x);

	// 隣接領域(edgeのさらに1つ外側)
	double *yup_surr = (double *)malloc(sizeof(double) * localGridSize_x);
	double *ydown_surr = (double *)malloc(sizeof(double) * localGridSize_x);

	int i = 0;
	for (i = 0; i < NITER; i++)
	{
		// 担当領域の値を計算
		for (int localx = 1; localx <= localGridSize_x; localx++)
		{
			for (int localy = 1; localy <= localGridSize_y; localy++)
			{
				u_new[localx][localy] = 0.25f * (u_old[localx - 1][localy - 1] + u_old[localx - 1][localy + 1] + u_old[localx + 1][localy - 1] + u_old[localx + 1][localy + 1]);
			}
		}

		// 周辺領域を同期
		for (int localx = 1; localx <= localGridSize_x; localx++)
		{
			yup_edge[localx - 1] = u_new[localx][gridSize_y];
			ydown_edge[localx - 1] = u_new[localx][1];
		}
		MPI_Irecv(yup_surr, localGridSize_x, MPI_DOUBLE, yup, 0, comm2d, &req_yup);
		MPI_Irecv(ydown_surr, localGridSize_x, MPI_DOUBLE, ydown, 0, comm2d, &req_ydown);
		MPI_Irecv(&(u_new[localGridSize_x + 1][1]), localGridSize_y, MPI_DOUBLE, xup, 0, comm2d, &req_xup);
		MPI_Irecv(&(u_new[0][1]), localGridSize_y, MPI_DOUBLE, xdown, 0, comm2d, &req_xdown);

		MPI_Send(yup_edge, localGridSize_x, MPI_DOUBLE, yup, 0, comm2d);
		MPI_Send(ydown_edge, localGridSize_x, MPI_DOUBLE, ydown, 0, comm2d);
		MPI_Send(&(u_new[localGridSize_x][1]), localGridSize_y, MPI_DOUBLE, xup, 0, comm2d);
		MPI_Send(&(u_new[1][1]), localGridSize_y, MPI_DOUBLE, xdown, 0, comm2d);
		for (int localx = 1; localx <= localGridSize_x; localx++)
		{
			u_new[localx][localGridSize_y + 1] = yup_surr[localx - 1];
			u_new[localx][0] = ydown_surr[localx - 1];
		}

		double **temp = u_old;
		u_old = u_new;
		u_new = temp;
	}

	// 総和を求める
	double localsum = 0;
	for (int localx = 1; localx <= localGridSize_x; localx++)
	{
		for (int localy = 1; localy <= localGridSize_y; localy++)
		{
			localsum += u_old[localx][localy];
		}
	}
	double sum = -100;
	MPI_Reduce(&localsum, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, comm2d);
	if (myrank == 0)
	{
		printf("sum = %g\n", sum);
	}
	MPI_Comm_free(comm2d);

	MPI_Barrier(comm2d);
	double end = MPI_Wtime();
	if (myrank == 0)
	{
		printf("time = %g\n", end - start);
	}
	MPI_Finalize();

	return (0);
}