/*
 * Laplace equation with explicit method
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <sys/stat.h>
/* square region */
#define PI 3.1415927

#define XSIZE 256
#define YSIZE 256
#define NITER 10000


//#define XSIZE 64
//#define YSIZE 64
//#define NITER 4

#define LOCALITER 4

#define TAG_FROM_XUP 800
#define TAG_FROM_XDOWN 801
#define TAG_FROM_YUP 810
#define TAG_FROM_YDOWN 811

#define TAG_FROM_XUP2 900
#define TAG_FROM_XDOWN2 901
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
	int xupEdge=LOCALITER+localGridSize_x-1;
	int yupEdge=LOCALITER+localGridSize_y-1;

	//debug
	//fprintf(stderr,"rank:%d pos(%d,%d) gridsize(%d,%d) xrange:%d-%d yrange:%d-%d localArySize_x:%d localArySize_y:%d localGridSize_x:%d localGridSize_y:%d\n", myrank, gridPos_x, gridPos_y, localGridSize_x, localGridSize_y,xdownEdge,xupEdge,ydownEdge,yupEdge,localArySize_x, localArySize_y,localGridSize_x,localGridSize_y);
	MPI_Barrier(comm2d);


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
			int worldy = localy - ydownEdge + ystart;
			//領域全体の外側に当たる部分の値は0にする。
			if (worldx < 0 || worldy < 0 || XSIZE - 1 < worldx || YSIZE - 1 < worldy)
			{
				u_old[localx][localy] = 0;
				u_new[localx][localy] = 0;
			}
			else
			{
				u_old[localx][localy] = sin((double)worldx / XSIZE * PI) + cos((double)worldy / YSIZE * PI);
				u_new[localx][localy] = u_old[localx][localy];
			}
		}
	}

	// laplaceを解く--------------------------------------------------------
	MPI_Barrier(MPI_COMM_WORLD);
	double start = MPI_Wtime();

	int xdown, xup, ydown, yup;
	int xdyd=-2, xdyu=-2, xuyd=-2, xuyu=-2;
	MPI_Request req_yuTxu,req_yuTxd,req_ydTxu,req_ydTxd;
	MPI_Status sta_yuTxu_R,sta_yuTxd_R,sta_ydTxu_R,sta_ydTxd_R;
	MPI_Status sta_yuTxu_S,sta_yuTxd_S,sta_ydTxu_S,sta_ydTxd_S;
	//上下左右で隣接するプロセスのランクを取得
	MPI_Cart_shift(comm2d, 0, 1, &xdown, &xup);
	MPI_Cart_shift(comm2d, 1, 1, &ydown, &yup);
	//右上　右下　左上　左下のプロセスのランクを他のプロセスから受け取る。
	MPI_Isend(&yup, 1, MPI_INT,xup, TAG_FROM_XUP, comm2d, &req_yuTxu);
	MPI_Isend(&yup, 1, MPI_INT,xdown, TAG_FROM_YDOWN, comm2d, &req_yuTxd);
	MPI_Isend(&ydown, 1, MPI_INT,xup, TAG_FROM_XUP2, comm2d, &req_ydTxu);
	MPI_Isend(&ydown, 1, MPI_INT,xdown, TAG_FROM_YDOWN2, comm2d, &req_ydTxd);
	MPI_Recv(&xdyu,1,MPI_INT,xdown,TAG_FROM_XUP,comm2d,&sta_yuTxu_R);
	MPI_Recv(&xuyu,1,MPI_INT,xup,TAG_FROM_YDOWN,comm2d,&sta_yuTxd_R);
	MPI_Recv(&xdyd,1,MPI_INT,xdown,TAG_FROM_XUP2,comm2d,&sta_ydTxu_R);
	MPI_Recv(&xuyd,1,MPI_INT,xup,TAG_FROM_YDOWN2,comm2d,&sta_ydTxd_R);
	MPI_Wait(&req_yuTxu, &sta_yuTxu_S);
	MPI_Wait(&req_yuTxd, &sta_yuTxd_S);
	MPI_Wait(&req_ydTxu, &sta_ydTxu_S);
	MPI_Wait(&req_ydTxd, &sta_ydTxd_S);

	//fprintf(stderr,"xdown:%d, xup:%d, ydown:%d, yup:%d, xdyd:%d, xdyu:%d, xuyd:%d, xuyu:%d\n",xdown, xup, ydown, yup,xdyd, xdyu, xuyd, xuyu);


	MPI_Request req_xup, req_xdown, req_yup, req_ydown;
	MPI_Request req_xdyd, req_xdyu, req_xuyd, req_xuyu;
	MPI_Status status_xup, status_xdown, status_yup, status_ydown;
	MPI_Status status_xdyd, status_xdyu, status_xuyd, status_xuyu;
	MPI_Status status_toxup, status_toxdown, status_toyup, status_toydown;
	MPI_Status status_toxdyd, status_toxdyu, status_toxuyd, status_toxuyu;

	//y方向の隣接領域を担当するプロセスとの通信に使うバッファを確保
	//x方向はデータを格納している領域を使って直接データを送受信する


	// 担当領域の一番端(edge)
	double **yup_edge_Buf = (double **)malloc(sizeof(double *) * LOCALITER);
	double **ydown_edge_Buf = (double **)malloc(sizeof(double *) * LOCALITER);
	// 隣接領域(edgeのさらに1つ外側)
	double **yup_surr_Buf = (double **)malloc(sizeof(double *) * LOCALITER);
	double **ydown_surr_Buf = (double **)malloc(sizeof(double *) * LOCALITER);

	int xDirHaloExch_ExchValNum=LOCALITER*localGridSize_y;
	int yDirHaloExch_ExchValNum=LOCALITER*localGridSize_x;
	// 隣接領域(edgeのさらに1つ外側)(Halo)の値を入れるバッファ
	double *xupHalo_ValBuf=malloc(sizeof(double)*xDirHaloExch_ExchValNum);
	double *xdownHalo_ValBuf=malloc(sizeof(double)*xDirHaloExch_ExchValNum);
	double *yupHalo_ValBuf=malloc(sizeof(double)*yDirHaloExch_ExchValNum);
	double *ydownHalo_ValBuf=malloc(sizeof(double)*yDirHaloExch_ExchValNum);

	int cornerHalo_BufSize=(LOCALITER-1+1)*(LOCALITER-1)/2;
	double *xdydHalo_ValBuf=malloc(sizeof(double)*cornerHalo_BufSize);
	double *xdyuHalo_ValBuf=malloc(sizeof(double)*cornerHalo_BufSize);
	double *xuydHalo_ValBuf=malloc(sizeof(double)*cornerHalo_BufSize);
	double *xuyuHalo_ValBuf=malloc(sizeof(double)*cornerHalo_BufSize);

	// 担当領域の一番端(edge)の値を入れるバッファ
	double *xupEdge_ValBuf=malloc(sizeof(double)*LOCALITER*localGridSize_y);
	double *xdownEdge_ValBuf=malloc(sizeof(double)*LOCALITER*localGridSize_y);
	double *yupEdge_ValBuf=malloc(sizeof(double)*LOCALITER*localGridSize_x);
	double *ydownEdge_ValBuf=malloc(sizeof(double)*LOCALITER*localGridSize_x);
	
	double *xdydEdge_ValBuf=malloc(sizeof(double)*cornerHalo_BufSize);
	double *xdyuEdge_ValBuf=malloc(sizeof(double)*cornerHalo_BufSize);
	double *xuydEdge_ValBuf=malloc(sizeof(double)*cornerHalo_BufSize);
	double *xuyuEdge_ValBuf=malloc(sizeof(double)*cornerHalo_BufSize);

	

	for(int i=0;i<LOCALITER;i++)
	{
		yup_edge_Buf[i]=(double *)malloc(sizeof(double)*localGridSize_x);
		ydown_edge_Buf[i]=(double *)malloc(sizeof(double)*localGridSize_x);
		yup_surr_Buf[i]=(double *)malloc(sizeof(double)*localGridSize_x);
		ydown_surr_Buf[i]=(double *)malloc(sizeof(double)*localGridSize_x);
	}
	//debug
	MPI_Barrier(comm2d);
	//fprintf(stderr,"rank:%d pos(%d,%d) gridsize(%d,%d) xrange:%d-%d yrange:%d-%d xup:%d xdown:%d yup:%d ydown:%d \n", myrank, gridPos_x, gridPos_y, localGridSize_x, localGridSize_y,xstart,xend,ystart,yend, xup, xdown, yup, ydown);


	for (int i = 0; i < NITER; i+=LOCALITER)
	{
		
		int restIter=NITER-i;
		int localIterNum=LOCALITER<restIter?LOCALITER:restIter;

		for(int calHaloDepth=localIterNum-1;calHaloDepth>=0;calHaloDepth--)
		{
			
			//このループで計算するべきHaloの深さがcalHaloDepth
			//今回計算する点のxの範囲を求める
			int lxFirst=(gridPos_x>0)?xdownEdge-calHaloDepth:xdownEdge;
			int lxLast=(gridPos_x<gridNum_x-1)?xupEdge+calHaloDepth:xupEdge;
			// 担当領域の値と、次のループの計算に使う周辺領域の値の計算を行う
			for (localx = lxFirst; localx <= lxLast; localx++)
			{
				//そのxにおけるy方向のHaloの大きさを求める
				int thisX_HaloYDepth=calHaloDepth;
				if(localx<xdownEdge)
				{
					thisX_HaloYDepth=lxFirst-localx;
				}
				else if(localx>xupEdge)
				{
					thisX_HaloYDepth=localx-lxLast;
				}
				
				//今回計算する点のyの範囲を求める
				//プロセス自体の位置が端だった場合、それより外側は計算しない。
				int lyFirst=(gridPos_y>0)?ydownEdge-thisX_HaloYDepth:ydownEdge;
				int lyLast=(gridPos_y<gridNum_y-1)?yupEdge+thisX_HaloYDepth:yupEdge;

				for (localy = lyFirst; localy <=lyLast; localy++)
				{
					//fprintf(stderr,"i:%d (%d,%d)(%d,%d) \n",i,gridPos_x,gridPos_y,localx,localy);
					u_new[localx][localy] = 0.25 * (u_old[localx - 1][localy] + u_old[localx + 1][localy] + u_old[localx][localy - 1] + u_old[localx][localy + 1]);
				}
			}
			double **temp=u_new;
			u_new=u_old;
			u_old=temp;
		}
		double **temp2=u_new;
		u_new=u_old;
		u_old=temp2;

		int diagEdgeCopy_StartBufIndex=0;
		//斜め、隣接している領域を担当するプロセスに送るデータをバッファに送る
		for (int exchEdgeDepth=1;exchEdgeDepth<=LOCALITER;exchEdgeDepth++)
		{
			int yEdgeCopy_StartBufIndex=(exchEdgeDepth-1)*localGridSize_x;
			//y方向の隣接領域を担当するプロセスに送るデータをバッファに入れる
			for(int localx=xdownEdge;localx<=xupEdge;localx++)
			{
				int yEdgeCopy_BufIndex=yEdgeCopy_StartBufIndex+localx-xdownEdge;
				yupEdge_ValBuf[yEdgeCopy_BufIndex]=u_new[localx][yupEdge-(exchEdgeDepth-1)];
				ydownEdge_ValBuf[yEdgeCopy_BufIndex]=u_new[localx][ydownEdge+(exchEdgeDepth-1)];
			}

			//x方向の隣接領域を担当するプロセスに送るデータをバッファに入れる
			int xEdgeCopy_StartBufIndex=(exchEdgeDepth-1)*localGridSize_y;
			memcpy(&(xupEdge_ValBuf[xEdgeCopy_StartBufIndex]), &(u_new[xupEdge-(exchEdgeDepth-1)][ydownEdge]), sizeof(double)*localGridSize_y);
			memcpy(&(xdownEdge_ValBuf[xEdgeCopy_StartBufIndex]), &(u_new[xdownEdge+(exchEdgeDepth-1)][ydownEdge]), sizeof(double)*localGridSize_y);
			
			
			//debug
			for(int localy=ydownEdge;localy<=yupEdge;localy++)
			{
				int xEdgeCopy_BufIndex=xEdgeCopy_StartBufIndex+localy-ydownEdge;
				if(xupEdge_ValBuf[xEdgeCopy_BufIndex]!=u_new[xupEdge-(exchEdgeDepth-1)][localy])
				{
					fprintf(stderr,"ffff\n");
				}
				if(xdownEdge_ValBuf[xEdgeCopy_BufIndex]!=u_new[xdownEdge+(exchEdgeDepth-1)][localy])
				{
					fprintf(stderr,"ssss\n");
				}	
			}
			
			

			//斜めにある領域を担当するプロセスに送るデータをバッファに入れる
			int diagEdgeCopy_CopyValNum=LOCALITER-exchEdgeDepth;
			memcpy(&(xdydEdge_ValBuf[diagEdgeCopy_StartBufIndex]),&(u_new[xdownEdge+(exchEdgeDepth-1)][ydownEdge]),sizeof(double)*diagEdgeCopy_CopyValNum);
			memcpy(&(xdyuEdge_ValBuf[diagEdgeCopy_StartBufIndex]),&(u_new[xdownEdge+(exchEdgeDepth-1)][yupEdge-(diagEdgeCopy_CopyValNum-1)]),sizeof(double)*diagEdgeCopy_CopyValNum);
			memcpy(&(xuydEdge_ValBuf[diagEdgeCopy_StartBufIndex]),&(u_new[xupEdge-(exchEdgeDepth-1)][ydownEdge]),sizeof(double)*diagEdgeCopy_CopyValNum);
			memcpy(&(xuyuEdge_ValBuf[diagEdgeCopy_StartBufIndex]),&(u_new[xupEdge-(exchEdgeDepth-1)][yupEdge-(diagEdgeCopy_CopyValNum-1)]),sizeof(double)*diagEdgeCopy_CopyValNum);
			diagEdgeCopy_StartBufIndex+=diagEdgeCopy_CopyValNum;
		}

		//2次元分割なので、上下左右のプロセスと同期する。
		//上下左右の端の領域を隣接プロセスに送信
		MPI_Request sendToXdown_Request,sendToXup_Request,sendToYdown_Request,sendToYup_Request;
		MPI_Isend(yupEdge_ValBuf, yDirHaloExch_ExchValNum, MPI_DOUBLE, yup, BASE_FROM_YDOWN+1, comm2d, &sendToYup_Request);
		MPI_Isend(ydownEdge_ValBuf, yDirHaloExch_ExchValNum, MPI_DOUBLE, ydown, BASE_FROM_YUP+1, comm2d, &sendToYdown_Request);
		MPI_Isend(xupEdge_ValBuf, xDirHaloExch_ExchValNum, MPI_DOUBLE, xup, BASE_FROM_XDOWN+1, comm2d, &sendToXup_Request);
		MPI_Isend(xdownEdge_ValBuf, xDirHaloExch_ExchValNum, MPI_DOUBLE, xdown, BASE_FROM_XUP+1, comm2d, &sendToXdown_Request);
		//4隅の端の領域を斜めにあるプロセスへ送信する
		MPI_Request sendToXdyd_Request,sendToXdyu_Request,sendToXuyd_Request,sendToXuyu_Request;
		MPI_Isend(xdydEdge_ValBuf,cornerHalo_BufSize,MPI_DOUBLE,xdyd,BASE_FROM_XUPYUP+1,comm2d,&sendToXdyd_Request);
		MPI_Isend(xdyuEdge_ValBuf,cornerHalo_BufSize,MPI_DOUBLE,xdyu,BASE_FROM_XUPYDOWN+1,comm2d,&sendToXdyu_Request);
		MPI_Isend(xuydEdge_ValBuf,cornerHalo_BufSize,MPI_DOUBLE,xuyd,BASE_FROM_XDOWNYUP+1,comm2d,&sendToXuyd_Request);
		MPI_Isend(xuyuEdge_ValBuf,cornerHalo_BufSize,MPI_DOUBLE,xuyu,BASE_FROM_XDOWNYDOWN+1,comm2d,&sendToXuyu_Request);

		
		//上下左右にあるプロセスの端の領域をそのプロセスから受信
		MPI_Status recvFromXdown_Status,  recvFromXup_Status, recvFromYdown_Status, recvFromYup_Status;
		if(ydown!=MPI_PROC_NULL){MPI_Recv(ydownHalo_ValBuf, yDirHaloExch_ExchValNum, MPI_DOUBLE, ydown, BASE_FROM_YDOWN+1, comm2d, &recvFromYdown_Status);}
		if(yup!=MPI_PROC_NULL){MPI_Recv(yupHalo_ValBuf, yDirHaloExch_ExchValNum, MPI_DOUBLE, yup, BASE_FROM_YUP+1, comm2d, &recvFromYup_Status);}
		if(xdown!=MPI_PROC_NULL){MPI_Recv(xdownHalo_ValBuf, xDirHaloExch_ExchValNum, MPI_DOUBLE, xdown, BASE_FROM_XDOWN+1, comm2d, &recvFromXdown_Status);}
		if(xup!=MPI_PROC_NULL){MPI_Recv(xupHalo_ValBuf, xDirHaloExch_ExchValNum, MPI_DOUBLE, xup, BASE_FROM_XUP+1, comm2d, &recvFromXup_Status);}
			
		//斜めにあるプロセスの端の領域をそのプロセスから受信
		MPI_Status recvFromXdyd_Status,  recvFromXdyu_Status, recvFromXuyd_Status, recvFromXuyu_Status;
		if(xuyu!=MPI_PROC_NULL){MPI_Recv(xuyuHalo_ValBuf,cornerHalo_BufSize,MPI_DOUBLE,xuyu,BASE_FROM_XUPYUP+1,comm2d,&recvFromXuyu_Status);}
		if(xuyd!=MPI_PROC_NULL){MPI_Recv(xuydHalo_ValBuf,cornerHalo_BufSize,MPI_DOUBLE,xuyd,BASE_FROM_XUPYDOWN+1,comm2d,&recvFromXuyd_Status);}
		if(xdyu!=MPI_PROC_NULL){MPI_Recv(xdyuHalo_ValBuf,cornerHalo_BufSize,MPI_DOUBLE,xdyu,BASE_FROM_XDOWNYUP+1,comm2d,&recvFromXdyu_Status);}
		if(xdyd!=MPI_PROC_NULL){MPI_Recv(xdydHalo_ValBuf,cornerHalo_BufSize,MPI_DOUBLE,xdyd,BASE_FROM_XDOWNYDOWN+1,comm2d,&recvFromXdyd_Status);}

		MPI_Status sendToXdown_Status,  sendToXup_Status, sendToYdown_Status, sendToYup_Status;
		MPI_Wait(&sendToXdown_Request, &sendToXdown_Status);
		MPI_Wait(&sendToXup_Request, &sendToXup_Status);
		MPI_Wait(&sendToYdown_Request, &sendToYdown_Status);
		MPI_Wait(&sendToYup_Request, &sendToYup_Status);

		MPI_Status sendToXdyd_Status,  sendToXdyu_Status, sendToXuyd_Status, sendToXuyu_Status;
		MPI_Wait(&sendToXdyd_Request, &sendToXdyd_Status);
		MPI_Wait(&sendToXdyu_Request, &sendToXdyu_Status);
		MPI_Wait(&sendToXuyd_Request, &sendToXuyd_Status);
		MPI_Wait(&sendToXuyu_Request, &sendToXuyu_Status);

		
		

		int diagHaloCopy_StartBufIndex=0;
		for(int EdgeOrHaloDepth=1;EdgeOrHaloDepth<=LOCALITER;EdgeOrHaloDepth++)
		{
			//x方向の隣接領域から送られてきたデータをu_newに代入していく。
			int xDirHaloCopy_ValBufStartIndex=(EdgeOrHaloDepth-1)*localGridSize_y;
			if(xdown!=MPI_PROC_NULL){memcpy(&(u_new[xdownEdge-EdgeOrHaloDepth][ydownEdge]),&(xdownHalo_ValBuf[xDirHaloCopy_ValBufStartIndex]),localGridSize_y*sizeof(double));}
			if(xup!=MPI_PROC_NULL){memcpy(&(u_new[xupEdge+EdgeOrHaloDepth][ydownEdge]),&(xupHalo_ValBuf[xDirHaloCopy_ValBufStartIndex]),localGridSize_y*sizeof(double));}

			//y方向の隣接領域から送られてきたデータをu_newに代入していく。
			int yDirHaloCopy_ValBufStartIndex=(EdgeOrHaloDepth-1)*localGridSize_x;
			for(int localx=xdownEdge;localx<=xupEdge;localx++)
			{
				int yDirHaloCopy_BufIndex=yDirHaloCopy_ValBufStartIndex+(localx-xdownEdge);
				u_new[localx][ydownEdge-EdgeOrHaloDepth]= (ydown != MPI_PROC_NULL) ? ydownHalo_ValBuf[yDirHaloCopy_BufIndex]:0;
				u_new[localx][yupEdge+EdgeOrHaloDepth]= (yup != MPI_PROC_NULL) ?yupHalo_ValBuf[yDirHaloCopy_BufIndex]:0;
			}
			//斜めの領域を担当するプロセスから送られてきたデータをu_newに代入していく。
			//斜めにある領域を担当するプロセスに送るデータをバッファに入れる
			int diagHaloCopy_CopyValNum=LOCALITER-EdgeOrHaloDepth;
			if(xdyd!=MPI_PROC_NULL){memcpy(&(u_new[xdownEdge-(EdgeOrHaloDepth)][ydownEdge-diagHaloCopy_CopyValNum]),&(xdydHalo_ValBuf[diagHaloCopy_StartBufIndex]),sizeof(double)*diagHaloCopy_CopyValNum);}
			if(xdyu!=MPI_PROC_NULL){memcpy(&(u_new[xdownEdge-(EdgeOrHaloDepth)][yupEdge+1]),&(xdyuHalo_ValBuf[diagHaloCopy_StartBufIndex]),sizeof(double)*diagHaloCopy_CopyValNum);}
			if(xuyd!=MPI_PROC_NULL){memcpy(&(u_new[xupEdge+(EdgeOrHaloDepth)][ydownEdge-diagHaloCopy_CopyValNum]),&(xuydHalo_ValBuf[diagHaloCopy_StartBufIndex]),sizeof(double)*diagHaloCopy_CopyValNum);}
			if(xuyu!=MPI_PROC_NULL){memcpy(&(u_new[xupEdge+(EdgeOrHaloDepth)][yupEdge+1]),&(xuyuHalo_ValBuf[diagHaloCopy_StartBufIndex]),sizeof(double)*diagHaloCopy_CopyValNum);}
			diagHaloCopy_StartBufIndex+=diagHaloCopy_CopyValNum;
		}
		//データコピーの代わりに、書き込む領域を入れ替える
		double **temp = u_old;
		u_old = u_new;
		u_new = temp;
	}

	double localsum = 0;
	for (localx = xdownEdge; localx <=xupEdge ; localx++)
	{
		for (localy = ydownEdge; localy <= yupEdge; localy++)
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

	//debug
	char dirName_str[256];
	sprintf(dirName_str, "log/Liter(%d,%d)_%d", gridNum_x,gridNum_y,NITER);
	mkdir(dirName_str, 0777);
	char fileName_str[256];
	sprintf(fileName_str, "log/Liter(%d,%d)_%d/(%d,%d)_datafile.txt", gridNum_x,gridNum_y,NITER,gridPos_x,gridPos_y);
	FILE * dataFile_pointer= fopen(fileName_str, "w") ;
	fprintf(dataFile_pointer,"   x\\y:");
	for(int localy=	0;localy<localArySize_y;localy++)
	{
		int worldy = localy - ydownEdge+ystart;
		if(ydownEdge<=localy&&localy<=yupEdge)
		{
			fprintf(dataFile_pointer,"%7d|",worldy);
		}
		else 
		{
			fprintf(dataFile_pointer,"(%5d)|",worldy);
		}
		
	}
	fprintf(dataFile_pointer,"\n");
	for(int localx=0;localx<localArySize_x;localx++)
	{
		int worldx=localx-xdownEdge+xstart;
		if(xdownEdge<=localx&&localx<=xupEdge)
		{
			fprintf(dataFile_pointer,"%6d:",worldx);
		}
		else
		{
			fprintf(dataFile_pointer,"(%+4d):",worldx);
		}

		for(int localy=	0;localy<localArySize_y;localy++)
		{
			fprintf(dataFile_pointer,"%+1.4lf|",u_old[localx][localy]);
		}

		if(xdownEdge<=localx&&localx<=xupEdge)
		{
			fprintf(dataFile_pointer,"%6d:",worldx);
		}
		else
		{
			fprintf(dataFile_pointer,"(%+4d):",worldx);
		}
		fprintf(dataFile_pointer,"\n");
	}

	fprintf(dataFile_pointer,"   x\\y:");
	for(int localy=	0;localy<localArySize_y;localy++)
	{
		int worldy = localy - ydownEdge+ystart;
		if(ydownEdge<=localy&&localy<=yupEdge)
		{
			fprintf(dataFile_pointer,"%7d|",worldy);
		}
		else 
		{
			fprintf(dataFile_pointer,"(%5d)|",worldy);
		}
		
	}
	fprintf(dataFile_pointer,"\n");
	fclose(dataFile_pointer);

	MPI_Finalize();
	return (0);
}