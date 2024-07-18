/*
 * Laplace equation with explicit method
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <sys/stat.h>
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
	fprintf(stderr,"rank:%d pos(%d,%d) gridsize(%d,%d) xrange:%d-%d yrange:%d-%d localArySize_x:%d localArySize_y:%d localGridSize_x:%d localGridSize_y:%d\n", myrank, gridPos_x, gridPos_y, localGridSize_x, localGridSize_y,xdownEdge,xupEdge,ydownEdge,yupEdge,localArySize_x, localArySize_y,localGridSize_x,localGridSize_y);
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
			//自分が担当する領域の一番端のindexがLocalIter
			//このループで計算するべきHaloの深さがcalHaloDepth
			//今回計算する点の範囲を求める

			
			
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
				//プロセス自体の位置が端だった場合、それより外側は計算しない。
				int lyFirst=(gridPos_y>0)?ydownEdge-thisX_HaloYDepth:ydownEdge;
				int lyLast=(gridPos_y<gridNum_y-1)?yupEdge+thisX_HaloYDepth:yupEdge;

				//debug
				//fprintf(stderr,"i:%d/%d (%d,%d) lx:%d ly%d~%d\n",i,calHaloDepth,gridPos_x,gridPos_y,localx,lyFirst,lyLast);

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

		// 自分の担当領域の境界部分を隣接プロセスと同期
		//y方向の隣接領域を担当するプロセスに送るデータをバッファに入れる
		for (int exchEdgeDepth=1;exchEdgeDepth<=LOCALITER;exchEdgeDepth++)
		{
			
			for(int localx=xdownEdge;localx<=xupEdge;localx++)
			{
				yup_edge_Buf[exchEdgeDepth-1][localx-xdownEdge] = u_new[localx][yupEdge-(exchEdgeDepth-1)];
				ydown_edge_Buf[exchEdgeDepth-1][localx-xdownEdge] = u_new[localx][ydownEdge+(exchEdgeDepth-1)];
			}
			
		}

		
		//深さexchHaloDepthのHaloを交換する
		for(int exchHaloDepth=1;exchHaloDepth<=LOCALITER;exchHaloDepth++)
		{	
			//2次元分割なので、上下左右のプロセスと同期する。
			//端の領域を隣接プロセスに送信
			MPI_Isend(yup_edge_Buf[exchHaloDepth-1], localGridSize_x, MPI_DOUBLE, yup, BASE_FROM_YDOWN+exchHaloDepth, comm2d, &req_yup);
			MPI_Isend(ydown_edge_Buf[exchHaloDepth-1], localGridSize_x, MPI_DOUBLE, ydown, BASE_FROM_YUP+exchHaloDepth, comm2d, &req_ydown);
			MPI_Isend(&(u_new[xupEdge-(exchHaloDepth-1)][ydownEdge]), localGridSize_y, MPI_DOUBLE, xup, BASE_FROM_XDOWN+exchHaloDepth, comm2d, &req_xup);
			MPI_Isend(&(u_new[xdownEdge+(exchHaloDepth-1)][ydownEdge]), localGridSize_y, MPI_DOUBLE, xdown, BASE_FROM_XUP+exchHaloDepth, comm2d, &req_xdown);
			//斜めにあるプロセスからは、x方向で深さexchHaloDepthのHaloを交換する
			//斜にある各プロセスと交換するHaloのサイズ
			int diagHaloNum=LOCALITER-(exchHaloDepth);
			MPI_Isend(&(u_new[xdownEdge+(exchHaloDepth-1)][0]),diagHaloNum,MPI_DOUBLE,xdyd,BASE_FROM_XUPYUP+exchHaloDepth,comm2d,&req_xdyd);
			MPI_Isend(&(u_new[xdownEdge+(exchHaloDepth-1)][localArySize_y-diagHaloNum]),diagHaloNum,MPI_DOUBLE,xdyu,BASE_FROM_XUPYDOWN+exchHaloDepth,comm2d,&req_xdyu);
			MPI_Isend(&(u_new[xupEdge-(exchHaloDepth-1)][0]),diagHaloNum,MPI_DOUBLE,xuyd,BASE_FROM_XDOWNYUP+exchHaloDepth,comm2d,&req_xuyd);
			MPI_Isend(&(u_new[xupEdge-(exchHaloDepth-1)][localArySize_y-diagHaloNum]),diagHaloNum,MPI_DOUBLE,xuyu,BASE_FROM_XDOWNYDOWN+exchHaloDepth,comm2d,&req_xuyu);

			
			//隣接プロセスの端の領域をほかプロセスから受信
			MPI_Status status_xup, status_xdown, status_yup, status_ydown;
			if(yup!=MPI_PROC_NULL){MPI_Recv(yup_surr_Buf[exchHaloDepth-1], localGridSize_x, MPI_DOUBLE, yup, BASE_FROM_YUP+exchHaloDepth, comm2d, &status_yup);}
			if(ydown!=MPI_PROC_NULL){MPI_Recv(ydown_surr_Buf[exchHaloDepth-1], localGridSize_x, MPI_DOUBLE, ydown, BASE_FROM_YDOWN+exchHaloDepth, comm2d, &status_ydown);}
			if(xup!=MPI_PROC_NULL){MPI_Recv(&(u_new[xupEdge+exchHaloDepth][ydownEdge]), localGridSize_y, MPI_DOUBLE, xup, BASE_FROM_XUP+exchHaloDepth, comm2d, &status_xup);}
			if(xdown!=MPI_PROC_NULL){MPI_Recv(&(u_new[xdownEdge-exchHaloDepth][ydownEdge]), localGridSize_y, MPI_DOUBLE, xdown, BASE_FROM_XDOWN+exchHaloDepth, comm2d, &status_xdown);}
			
			//斜めにあるプロセスの端の領域をほかプロセスから受信
			if(xdyd!=MPI_PROC_NULL){MPI_Recv(&(u_new[xdownEdge-exchHaloDepth][ydownEdge-diagHaloNum]),diagHaloNum,MPI_DOUBLE,xdyd,BASE_FROM_XDOWNYDOWN+exchHaloDepth,comm2d,&status_xdyd);}
			if(xdyu!=MPI_PROC_NULL){MPI_Recv(&(u_new[xdownEdge-exchHaloDepth][yupEdge+1]),diagHaloNum,MPI_DOUBLE,xdyu,BASE_FROM_XDOWNYUP+exchHaloDepth,comm2d,&status_xdyu);}
			if(xuyd!=MPI_PROC_NULL){MPI_Recv(&(u_new[xupEdge+exchHaloDepth][ydownEdge-diagHaloNum]),diagHaloNum,MPI_DOUBLE,xuyd,BASE_FROM_XUPYDOWN+exchHaloDepth,comm2d,&status_xuyd);}
			if(xuyu!=MPI_PROC_NULL){MPI_Recv(&(u_new[xupEdge+exchHaloDepth][yupEdge+1]),diagHaloNum,MPI_DOUBLE,xuyu,BASE_FROM_XUPYUP+exchHaloDepth,comm2d,&status_xuyu);}

			MPI_Wait(&req_xup, &status_toxup);
			MPI_Wait(&req_xdown, &status_toxdown);
			MPI_Wait(&req_yup, &status_toyup);
			MPI_Wait(&req_ydown, &status_toydown);

			MPI_Wait(&req_xdyd, &status_toxdyd);
			MPI_Wait(&req_xdyu, &status_toxdyu);
			MPI_Wait(&req_xuyd, &status_toxuyd);
			MPI_Wait(&req_xuyu, &status_toxuyu);




			int xupCount=-100,xdownCount=-100,yupCount=-100,ydownCount=-100;
			int xuyuCount=-100,xdyuCount=-100,xuydCount=-100,xdydCount=-100;
			int retxup= MPI_Get_count(&status_xup,MPI_DOUBLE,&xupCount);
			int retxdown= MPI_Get_count(&status_xdown,MPI_DOUBLE,&xdownCount);
			int retyup= MPI_Get_count(&status_yup,MPI_DOUBLE,&yupCount);
			int retydown= MPI_Get_count(&status_ydown,MPI_DOUBLE,&ydownCount);

			int retxuyu= MPI_Get_count(&status_xuyu,MPI_DOUBLE,&xuyuCount);
			int retxdyu= MPI_Get_count(&status_xdyu,MPI_DOUBLE,&xdyuCount);
			int retxuyd= MPI_Get_count(&status_xuyd,MPI_DOUBLE,&xuydCount);
			int retxdyd= MPI_Get_count(&status_xdyd,MPI_DOUBLE,&xdydCount);
			
			if(xuyuCount!=diagHaloNum&&xuyu!=MPI_PROC_NULL)
			{
				fprintf(stderr,"bufWrn i:%d r:%d (%d+1,%d+1) as:%d ac:%d ERR:%d ERR:%d\n",i,myrank,gridPos_x,gridPos_y,diagHaloNum,xuyuCount,status_xup.MPI_ERROR,retxuyu);
			}
			if(xdyuCount!=diagHaloNum&&xdyu!=MPI_PROC_NULL)
			{
				fprintf(stderr,"bufWrn i:%d r:%d (%d-1,%d+1) as:%d ac:%d ERR:%d ERR:%d\n",i,myrank,gridPos_x,gridPos_y,diagHaloNum,xdyuCount,status_xdown.MPI_ERROR,retxdyu);
			}
			if(xuydCount!=diagHaloNum&&xuyd!=MPI_PROC_NULL)
			{
				fprintf(stderr,"bufWrn i:%d r:%d (%d+1,%d-1) as:%d ac:%d ERR:%d ERR:%d\n",i,myrank,gridPos_x,gridPos_y,diagHaloNum,xuydCount,status_yup.MPI_ERROR,retxuyd);
			}
			if(xdydCount!=diagHaloNum&&xdyd!=MPI_PROC_NULL)
			{
				fprintf(stderr,"bufWrn i:%d r:%d (%d-1,%d-1) as:%d ac:%d ERR:%d ERR:%d\n",i,myrank,gridPos_x,gridPos_y,diagHaloNum,xdydCount,status_ydown.MPI_ERROR,retxdyd);
			}
			

			if(xupCount!=localGridSize_y&&xup!=MPI_PROC_NULL)
			{
				fprintf(stderr,"bufWrn i:%d r:%d (%d+1,%d) as:%d ac:%d ERR:%d ERR:%d\n",i,myrank,gridPos_x,gridPos_y,localGridSize_y,xupCount,status_xup.MPI_ERROR,retxup);
			}
			if(xdownCount!=localGridSize_y&&xdown!=MPI_PROC_NULL)
			{
				fprintf(stderr,"bufWrn i:%d r:%d (%d-1,%d) as:%d ac:%d ERR:%d ERR:%d\n",i,myrank,gridPos_x,gridPos_y,localGridSize_y,xdownCount,status_xdown.MPI_ERROR,retxdown);
			}
			if(yupCount!=localGridSize_x&&yup!=MPI_PROC_NULL)
			{
				fprintf(stderr,"bufWrn i:%d r:%d (%d,%d+1) as:%d ac:%d ERR:%d ERR:%d\n",i,myrank,gridPos_x,gridPos_y,localGridSize_x,yupCount,status_yup.MPI_ERROR,retyup);
			}
			if(ydownCount!=localGridSize_x&&ydown!=MPI_PROC_NULL)
			{
				fprintf(stderr,"bufWrn i:%d r:%d (%d,%d-1) as:%d ac:%d ERR:%d ERR:%d\n",i,myrank,gridPos_x,gridPos_y,localGridSize_x,ydownCount,status_ydown.MPI_ERROR,retydown);
			}
		}

		//y方向の隣接領域から送られてきたデータをu_newに代入していく。
		for(int EdgeOrHaloDepth=1;EdgeOrHaloDepth<=LOCALITER;EdgeOrHaloDepth++)
		{
			for(int localx=xdownEdge;localx<=xupEdge;localx++)
			{
				u_new[localx][ydownEdge-EdgeOrHaloDepth]= (ydown != MPI_PROC_NULL) ? ydown_surr_Buf[EdgeOrHaloDepth-1][localx-xdownEdge]:0;
				u_new[localx][yupEdge+EdgeOrHaloDepth]= (yup != MPI_PROC_NULL) ?yup_surr_Buf[EdgeOrHaloDepth-1][localx-xdownEdge]:0;
			}
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
	fprintf(dataFile_pointer,"x\\y:");
	for(int localy=	ydownEdge;localy<=yupEdge;localy++)
	{
		int worldy = localy - ydownEdge+ystart;
		fprintf(dataFile_pointer,"%7d|",worldy);
	}
	fprintf(dataFile_pointer,"\n");
	for(int localx=xdownEdge;localx<=xupEdge;localx++)
	{
		int worldx=localx-xdownEdge+xstart;
		fprintf(dataFile_pointer,"%3d:",worldx);
		for(int localy=	ydownEdge;localy<=yupEdge;localy++)
		{
			fprintf(dataFile_pointer,"%+1.4lf|",u_old[localx][localy]);
		}
		fprintf(dataFile_pointer,"\n");
	}
	fclose(dataFile_pointer);

	MPI_Finalize();
	return (0);
}