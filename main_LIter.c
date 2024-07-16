/*
 * Laplace equation with explicit method
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
/* square region */
#define PI 3.1415927

//#define XSIZE 256
//#define YSIZE 256
//#define NITER 10000


#define XSIZE 64
#define YSIZE 64
#define NITER 10

#define LOCALITER 2

#define TAG_FROM_XUP 800
#define TAG_FROM_XDOWN 801
#define TAG_FROM_YUP 810
#define TAG_FROM_YDOWN 811

#define TAG_FROM_XUP2 900
#define TAG_FROM_XDOWN2 9601
#define TAG_FROM_YUP2 910
#define TAG_FROM_YDOWN2 911

#define BASE_FROM_XUP 10000
#define BASE_FROM_XDOWN 20000
#define BASE_FROM_YUP 30000
#define BASE_FROM_YDOWN 40000

#define BASE_FROM_XUPYUP 130000
#define BASE_FROM_XDOWNYUP 230000
#define BASE_FROM_XUPYDOWN 140000
#define BASE_FROM_XDOWNYDOWN 240000

#ifndef FALSE
#define FALSE 0
#endif

#define DIM 2

int main(int argc, char *argv[])
{
	int localx,localy,i;

	int namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	int myrank, numprocs;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Get_processor_name(processor_name, &namelen);
	//fprintf(stderr, "Process %d on %s\n", myrank, processor_name);

	// 各パラメーターの決定--------------------------------------
	//  x,y軸方向に幾つ分割されているか
	int gridNum_x;
	int gridNum_y;

	// 1グリッドが担当する範囲の大きさ
	int gridSize_x;
	int gridSize_y;

	//gridNum_x*gridNum_y=numprocsを満たす組み合わせの中で、gridNum_xとgridNum_yの差が最も小さくなるような組み合わせを探す。
	gridNum_x = (int)(sqrtf(numprocs) + 1);
	for (gridNum_x = (int)(sqrtf(numprocs) + 1); gridNum_x >= 1; gridNum_x--)
	{
		if (numprocs % gridNum_x == 0)
		{
			gridNum_y = numprocs / gridNum_x;
			break;
		}
	}

	//実際の各gridの辺の長さ
	int *gridSizes_x=malloc(sizeof(int)*gridNum_x);
	int *gridSizes_y=malloc(sizeof(int)*gridNum_y);
	//各グリッドが始まる座標
	int *gridStart_x=malloc(sizeof(int)*gridNum_x);
	int *gridStart_y=malloc(sizeof(int)*gridNum_y);

	//グリッドの分割数から、グリッドの1辺の長さを求める
	gridSize_x = (XSIZE) / gridNum_x;
	gridSize_y = (YSIZE) / gridNum_y;
	//グリッドを敷き詰めても残ってしまう部分の辺の長さを求める。
	int restSize_x=XSIZE%gridNum_x;
	int restSize_y=YSIZE%gridNum_y;
;
	//余った辺の長さをそれぞれのグリッドに配って埋める。
	for(int gx=0;gx<gridNum_x;gx++)
	{
		if(gx<restSize_x)
		{
			gridSizes_x[gx]=gridSize_x+1;
			gridStart_x[gx]=gridSize_x*gx+gx;
		}
		else
		{
			gridSizes_x[gx]=gridSize_x;
			gridStart_x[gx]=gridSize_x*gx+restSize_x;
		}
	}
	for(int gy=0;gy<gridNum_y;gy++)
	{
		if(gy<restSize_y)
		{
			gridSizes_y[gy]=gridSize_y+1;
			gridStart_y[gy]=gridSize_y*gy+gy;
		}
		else
		{
			gridSizes_y[gy]=gridSize_y;
			gridStart_y[gy]=gridSize_y*gy+restSize_y;
		}
	}



	//2次元のキューブのトポロジーを持つコミュニケータを作成
	MPI_Comm comm2d;
	int periods[DIM] = {FALSE, FALSE};
	int dims[DIM] = {gridNum_x, gridNum_y};
	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, FALSE, &comm2d);


	//自プロセスがトポロジーのどの位置に割り当てられているかを取得
	int coords[DIM];
	int maxDim;
	MPI_Cart_coords(comm2d, myrank, 2, coords);
	int gridPos_x = coords[0];
	int gridPos_y = coords[1];

	// x,yそれぞれのgridSizeの総和がXSIZE,YSIZEよりも大きくなる場合があるので、
	//はみ出した分、端にあるグリッドの実際のgridSizeを小さくする。

	int xstart, xend, ystart, yend;
	xstart=gridStart_x[gridPos_x];
	xend=xstart+gridSizes_x[gridPos_x];
	ystart=gridStart_y[gridPos_y];
	yend=ystart+gridSizes_y[gridPos_y];
	
	int localGridSize_x = gridSizes_x[gridPos_x];
	int localGridSize_y = gridSizes_y[gridPos_y];
	//袖領域の値を格納する部分も含んでいる配列のサイズ
	int localArySize_x=localGridSize_x + LOCALITER*2;
	int localArySize_y=localGridSize_y + LOCALITER*2;
	//配列内の、自分の担当領域のx/y方向の端(edge)のx/y座標
	int xdownEdge=LOCALITER;
	int ydownEdge=LOCALITER;
	int xupEdge=LOCALITER+localArySize_x-1;
	int yupEdge=LOCALITER+localArySize_y-1;

	// uoldとunewのメモリ領域確保 + u_oldの初期化 ----------------------------
	//  u_oldとu_newは自分の担当領域+その隣接領域の値を保持する。
	//  担当するxの範囲は gridSize_x*gridPos_?-1 <= x <= gridSize_x*(gridPos_x+1) yの範囲も同様
	double **u_old = (double **)malloc(sizeof(double *) * (localArySize_x));
	double **u_new = (double **)malloc(sizeof(double *) * (localArySize_x));

	//localx,localyはu_old,u_newのインデックス
	for (localx = 0; localx < localArySize_x; localx++)
	{

		u_old[localx] = (double *)malloc(sizeof(double) * (localArySize_y));
		u_new[localx] = (double *)malloc(sizeof(double) * (localArySize_y));

		//worldx,worldyは分割される前の領域全体のインデックス
		int worldx = localx-xdownEdge+xstart;

		for (localy = 0; localy < localArySize_y; localy++)
		{
			int worldy = localy - ydown_edge + ystart;
			//領域全体の外側に当たる部分の値は0にする。
			if (worldx < 0 || worldy < 0 || XSIZE - 1 < worldx || YSIZE - 1 < worldy)
			{
				u_old[localx][localy] = 0;
				u_new[localx][localy] = 0;
			}
			else
			{
				u_old[localx][localy] = sin((double)worldx / XSIZE * PI) + cos((double)worldy / YSIZE * PI);
			}
		}
	}

	// laplaceを解く--------------------------------------------------------
	MPI_Barrier(MPI_COMM_WORLD);
	double start = MPI_Wtime();

	int xdown, xup, ydown, yup;
	int xdyd, xdyu, xuyd, xuyu;
	MPI_Request req_yuTxu,req_yuTxd,req_ydTxu,req_ydTxd;
	MPI_Status sta_yuTxu,sta_yuTxd_R,sta_ydTxu_R,sta_ydTxd_R;
	//上下左右で隣接するプロセスのランクを取得
	MPI_Cart_shift(comm2d, 0, 1, &xdown, &xup);
	MPI_Cart_shift(comm2d, 1, 1, &ydown, &yup);
	//右上　右下　左上　左下のプロセスのランクを他のプロセスから受け取る。
	MPI_Isend(&yup, 1, MPI_INT,xup, TAG_FROM_XUP, comm2d, &req_yuTxu);
	MPI_Isend(&yup, 1, MPI_INT,xdown, TAG_FROM_YDOWN, comm2d, &req_yuTxd);
	MPI_Isend(&ydown, 1, MPI_INT,xup, TAG_FROM_XUP2, comm2d, &req_ydTxu);
	MPI_Isend(&ydown, 1, MPI_INT,xdown, TAG_FROM_YDOWN2, comm2d, &req_ydTxd);
	MPI_Recv(&xdyu,1,MPI_INT,xdown,TAG_FROM_XUP,comm2d,sta_yuTxu);
	MPI_Recv(&xdyu,1,MPI_INT,xdown,TAG_FROM_XUP,comm2d,sta_yuTxu);
	MPI_Recv(&xdyu,1,MPI_INT,xdown,TAG_FROM_XUP,comm2d,sta_yuTxu);
	MPI_Recv(&xdyu,1,MPI_INT,xdown,TAG_FROM_XUP,comm2d,sta_yuTxu);


	MPI_Request req_xup, req_xdown, req_yup, req_ydown;
	MPI_Request req_xdyd, req_xdyu, req_xuyd, req_xuyu;
	MPI_Status status_xup, status_xdown, status_yup, status_ydown;
	MPI_Status status_xdyd, status_xdyu, status_xuyd, status_xuyu;
	MPI_Status status_toxup, status_toxdown, status_toyup, status_toydown;
	MPI_Status status_toxdyd, status_toxdyu, status_toxuyd, status_toxuyu;

	//y方向の隣接領域を担当するプロセスとの通信に使うバッファを確保
	//x方向はデータを格納している領域を使って直接データを送受信する

	// 担当領域の一番端(edge)
	double **yup_edge_Buf = (double **)malloc(sizeof(double*) * LOCALITER);
	double **ydown_edge_Buf = (double **)malloc(sizeof(double*) * LOCALITER);
	// 隣接領域(edgeのさらに1つ外側)
	double **yup_surr_Buf = (double **)malloc(sizeof(double*) * LOCALITER);
	double **ydown_surr_Buf = (double **)malloc(sizeof(double*) * LOCALITER);
	
	for(int i=0;i=LOCALITER;i++)
	{
		yup_edge_Buf[i]=(double *)malloc(sizeof(double)*localGridSize_x);
		ydown_edge_Buf[i]=(double *)malloc(sizeof(double)*localGridSize_x);
		yup_surr_Buf[i]=(double *)malloc(sizeof(double)*localGridSize_x);
		ydown_surr_Buf[i]=(double *)malloc(sizeof(double)*localGridSize_x);
	}

	for (i = 0; i < NITER; i++)
	{
		int restIter=NITER-i;
		int localIterNum=LOCALITER<restIter?LOCALITER:restIter;

		for(int calHaloDepth=localIterNum-1;calHaloDepth>=0;calHaloDepth--)
		{
			//自分が担当する領域の一番端のindexがLocalIter
			//このループで計算するべきHaloの深さがcalHaloDepth
			//今回計算する点の範囲を求める
			int lxFirst=xdownEdge-calHaloDepth;
			int lxLast=xupEdge+calHaloDepth;
			int lyFirst=ydown_edge-calHaloDepth;
			int lyLast=yup_edge+calHaloDepth;

			// 担当領域の値と、次のループの計算に使う周辺領域の値の計算を行う
			for (localx = lxFirst; localx <= lxLast; localx++)
			{
				for (localy = lyFirst; localy <=lyLast; localy++)
				{
					u_new[localx][localy] = 0.25 * (u_old[localx - 1][localy] + u_old[localx + 1][localy] + u_old[localx][localy - 1] + u_old[localx][localy + 1]);
				}
			}
		}

		// 自分の担当領域の境界部分を隣接プロセスと同期
		//y方向の隣接領域を担当するプロセスに送るデータをバッファに入れる
		for (int exchEdgeDepth=1;exchEdgeDepth<=LOCALITER;exchEdgeDepth++)
		{
			for(int localx=lxFirst;localx<lxLast;localx++)
			{
				yup_edge_Buf[exchEdgeDepth-1][localx-lxFirst] = u_new[localx][yupEdge-(exchEdgeDepth-1)];
				ydown_edge_Buf[exchEdgeDepth-1][localx-lxFirst] = u_new[localx][ydownEdge+(exchEdgeDepth-1)];
			}
			
		}

		
		//深さexchHaloDepthのHaloを交換する
		for(int exchHaloDepth=1;exchHaloDepth<=LOCALITER;exchHaloDepth++)
		{	
			
			//2次元分割なので、上下左右のプロセスと同期する。
			//端の領域を隣接プロセスに送信
			MPI_Isend(yup_edge_Buf[exchHaloDepth-1], localGridSize_x, MPI_DOUBLE, yup, BASE_FROM_YDOWN+exchHaloDepth, comm2d, &req_yup);
			MPI_Isend(ydown_edge_Buf[exchHaloDepth-1], localGridSize_x, MPI_DOUBLE, ydown, BASE_FROM_YUP+exchHaloDepth, comm2d, &req_ydown);
			MPI_Isend(&(u_new[localGridSize_x][1]), localGridSize_y, MPI_DOUBLE, xup, BASE_FROM_XDOWN+exchHaloDepth, comm2d, &req_xup);
			MPI_Isend(&(u_new[1][1]), localGridSize_y, MPI_DOUBLE, xdown, BASE_FROM_XUP+exchHaloDepth, comm2d, &req_xdown);
			//斜めにあるプロセスからは、x方向で深さexchHaloDepthのHaloを交換する
			//斜にある各プロセスと交換するHaloのサイズ
			int diagHaloNum=LOCALITER-(exchHaloDepth-1);
			MPI_Isend(&(u_new[xdownEdge+exchHaloDepth-1][0]),diagHaloNum,MPI_DOUBLE,xdyd,BASE_FROM_XUPYUP+exchHaloDepth,comm2d,&req_xdyd);
			MPI_Isend(&(u_new[xdownEdge+exchHaloDepth-1][localArySize_y-diagHaloNum]),diagHaloNum,MPI_DOUBLE,xdyu,BASE_FROM_XUPYDOWN+exchHaloDepth,comm2d,&req_xdyu);
			MPI_Isend(&(u_new[xupEdge-exchHaloDepth+1][0]),diagHaloNum,MPI_DOUBLE,xdyd,BASE_FROM_XUPYUP+exchHaloDepth,comm2d,&req_xdyd);
			MPI_Isend(&(u_new[xupEdge-exchHaloDepth+1][localArySize_y-exchHaloDepth]),diagHaloNum,MPI_DOUBLE,xdyu,BASE_FROM_XUPYDOWN+exchHaloDepth,comm2d,&req_xdyu);

			//隣接プロセスの端の領域をほかプロセスから受信
			MPI_Status status_xup, status_xdown, status_yup, status_ydown;
			if(yup!=MPI_PROC_NULL){MPI_Recv(yup_surr_Buf[exchHaloDepth-1], localGridSize_x, MPI_DOUBLE, yup, BASE_FROM_YUP+exchHaloDepth, comm2d, &status_yup);}
			if(ydown!=MPI_PROC_NULL){MPI_Recv(ydown_surr_Buf[exchHaloDepth-1], localGridSize_x, MPI_DOUBLE, ydown, BASE_FROM_YDOWN+exchHaloDepth, comm2d, &status_ydown);}
			if(xup!=MPI_PROC_NULL){MPI_Recv(&(u_new[localGridSize_x + 1][1]), localGridSize_y, MPI_DOUBLE, xup, BASE_FROM_XUP+exchHaloDepth, comm2d, &status_xup);}
			if(xdown!=MPI_PROC_NULL){MPI_Recv(&(u_new[0][1]), localGridSize_y, MPI_DOUBLE, xdown, BASE_FROM_XDOWN+exchHaloDepth, comm2d, &status_xdown);}
			//斜めにあるプロセスの端の領域をほかプロセスから受信
			if(xdyd!=MPI_PROC_NULL){ MPI_Recv(&(u_new[xdownEdge-exchHaloDepth][ydownEdge-diagHaloNum]),diagHaloNum,MPI_DOUBLE,xdyd,BASE_FROM_XDOWNYDOWN+exchHaloDepth,comm2d,status_xdyd);}
			if(xdyu!=MPI_PROC_NULL){MPI_Recv(&(u_new[xdownEdge-exchHaloDepth][yupEdge+1]),diagHaloNum,MPI_DOUBLE,xdyu,BASE_FROM_XDOWNYUP+exchHaloDepth,comm2dstatus_xdyu);}
			if(xuyd!=MPI_PROC_NULL){MPI_Recv(&(u_new[xupEdge+exchHaloDepth][ydownEdge-diagHaloNum]),diagHaloNum,MPI_DOUBLE,xuyd,BASE_FROM_XUPYDOWN+exchHaloDepth,comm2d,status_xuyd);}
			if(xuyu!=MPI_PROC_NULL){MPI_Recv(&(u_new[xupEdge+exchHaloDepth][yupEdge+1]),diagHaloNum,MPI_DOUBLE,xuyu,BASE_FROM_XUPYUP+exchHaloDepth,comm2d,status_xuyu);}

			MPI_Wait(&req_xup, &status_toxup);
			MPI_Wait(&req_xdown, &status_toxdown);
			MPI_Wait(&req_yup, &status_toyup);
			MPI_Wait(&req_ydown, &status_toydown);

			MPI_Wait(&req_xdyd, &status_toxdyd);
			MPI_Wait(&req_xdyu, &status_toxdyu);
			MPI_Wait(&req_xuyd, &status_toxuyd);
			MPI_Wait(&req_xuyu, &status_toxuyu);
		}

		//y方向の隣接領域から送られてきたデータをu_newに代入していく。
		for(int EdgeOrHaloDepth=1;EdgeOrHaloDepth<=LOCALITER;EdgeOrHaloDepth++)
		{
			for(int localx=lxFirst:localx<=lxLast;localx++)
			{
				u_new[localx][lyFirst-EdgeOrHaloDepth]=ydown_surr_Buf[EdgeOrHaloDepth-1][localx-lxFirst];
				u_new[localx][lyEnd+EdgeOrHaloDepth]=yup_surr_Buf[EdgeOrHaloDepth-1][localx-lxFirst];
			}
		}

		//データコピーの代わりに、書き込む領域を入れ替える
		double **temp = u_old;
		u_old = u_new;
		u_new = temp;
	}

	double localsum = 0;
	for (localx = 1; localx <= localGridSize_x; localx++)
	{
		for (localy = 1; localy <= localGridSize_y; localy++)
		{
			localsum += u_new[localx][localy]-u_old[localx][localy];
		}
	}
	
	double sum = -100;
	MPI_Reduce(&localsum, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, comm2d);
	if (myrank == 0)
	{
		fprintf(stderr,"sum = %lf\n", sum);
	}

	MPI_Comm_free(&comm2d);

	MPI_Barrier(MPI_COMM_WORLD);
	double end = MPI_Wtime();
	if (myrank == 0)
	{
		fprintf(stderr,"time = %g\n", end - start);
		printf("%g", end - start);
	}

	MPI_Finalize();
	return (0);
}