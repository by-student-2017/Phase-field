//多結晶粒組織形成のプログラム
//粒番号は1-nm
//粒番号nmを液相

#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

#define DRND(x) ((double)(x)/RAND_MAX*rand())//乱数の設定

//#define ND 100					//差分計算における計算領域一辺の分割数
//#define N 5							//考慮する結晶方位の数＋１

	//int nd=ND, ndm=ND-1;	//計算領域の一辺の差分分割数(差分ブロック数), ND-1を定義
	//int nm=N-1, nmm=N-2;	//考慮する結晶方位の数、N-2（考慮する結晶方位の数－１）を定義
	double PI=3.141592;		//pi
	double RR=8.3145;		//ガス定数
	double time1;			//計算カウント数
	//double ph[N][ND][ND], ph2[N][ND][ND];	//フェーズフィールド、フェーズフィールド補助配列
	//double aij[N][N];//勾配エネルギー係数
	//double wij[N][N];//ペナルティー項の係数
	//double tij[N][N];//粒界の易動度
	//double eij[N][N];//粒界移動の駆動力
	//int m00h[N][ND][ND];//位置(i,j)およびその周囲(i±1,j±1)において、pが０ではない方位の番号
	//int n00h[ND][ND];//位置(i,j)およびその周囲(i±1,j±1)において、pが０ではない方位の個数
	int Nstep, iout;

	void ini000(double *ph, int ND, int N);	//初期場の設定サブルーチン
	void datsave(double *ph, int ND, int N);	//データ保存サブルーチン
	void datsave_paraview(double *ph, int ND, int N);	//データ保存サブルーチン

//******* メインプログラム ******************************************
int main(void)
{
	int ND, N;
	int nd, ndm, nm, nmm;

	int i, j, k, l, ii, jj, kk, ll, it;	//整数
	int ip, im, jp, jm;									//整数
	int n1, n2, n3;											//整数
	int n00;		//位置(i,j)およびその周囲(i±1,j±1)において、pが０ではない方位の個数
	//int n000;		//位置(i,j)において、pが０ではない方位の個数（n00>=n000）
	double time1max;		//計算カウント数の最大値（計算終了カウント）
	double delt, L, b1;		//時間きざみ、計算領域１辺の長さ、差分ブロック一辺の長さ
	double M1, M0;			//粒界の易動度
	double W1, W0;			//ペナルティー項の係数
	double K1, K0;			//勾配エネルギー係数
	double E1, E0;			//粒界移動の駆動力
	double temp;			//温度
	double sum1, sum2, sum3;//各種の和の作業変数
	double pddtt;			//フェーズフィールドの時間変化率

	double gamma, gamma0;	//粒界エネルギ密度
	double delta;			//粒界幅（差分ブロック数にて表現）
	double amobi;			//粒界の易動度
	double vm0;				//モル体積

//****** 計算条件および物質定数の設定 ****************************************
	printf("---------------------------------\n");
	printf("read parameters from parameters.txt\n");
	FILE *fp;
	char name[40], comment[172];
	float param;
	float data[40];
	i = 0;
	fp = fopen("parameters.txt","r");
	while(fscanf(fp,"%s %e %[^\n] ",name, &param, comment) != EOF)
	{
		data[i] = param;
		printf("%d, %s %e \n",i,name,data[i]);
		i = i + 1;
	}
	fclose(fp);
	//
	ND      = int(data[0]);
	N       = int(data[1]);
	delt    = data[2];
	temp    = data[3];
	L       = data[4];
	vm0     = data[5];
	gamma0  = data[6];
	delta   = data[7];
	K0      = data[8];
	W0      = data[9];
	amobi   = data[10];
	M0      = data[11];
	E0      = data[12];
	time1max= int(data[13]);
	Nstep   = int(data[14]);
	printf("---------------------------------\n");
	//
	nd=ND;
	ndm=ND-1;
	nm=N-1;
	nmm=N-2;
	//
	double *ph   = (double *)malloc(sizeof(double)*( N*ND*ND + ND*ND + ND ));	//フェーズフィールド
	double *ph2  = (double *)malloc(sizeof(double)*( N*ND*ND + ND*ND + ND ));	//フェーズフィールド補助配列
	double *aij  = (double *)malloc(sizeof(double)*( ND*ND + ND ));	//勾配エネルギー係数
	double *wij  = (double *)malloc(sizeof(double)*( ND*ND + ND ));	//ペナルティー項の係数
	double *tij  = (double *)malloc(sizeof(double)*( ND*ND + ND ));	//粒界の易動度
	double *eij  = (double *)malloc(sizeof(double)*( ND*ND + ND ));	//粒界移動の駆動力
	double *m00h = (double *)malloc(sizeof(double)*( N*ND*ND + ND*ND + ND ));	//位置(i,j)およびその周囲(i±1,j±1)において、pが0ではない方位の番号
	double *n00h = (double *)malloc(sizeof(double)*( ND*ND + ND ));	//位置(i,j)およびその周囲(i±1,j±1)において、pが0ではない方位の個数
	//
	//printf("delt(5.0)=  "); scanf(" %lf",&delt);//時間刻みの入力	//	delt=5.0;

	//temp=1000.0; 					//温度(K)
	//L=2000.0;					//計算領域の一辺の長さ(nm)
	b1=L/(double)ND*1.0e-9; 	//差分プロック１辺の長さ(m)

	//vm0=7.0e-6;					/モル体積
	//gamma=0.5*vm0/RR/temp/b1;	//粒界エネルギ密度（0.5J/m^2）を無次元化
	gamma=gamma0*vm0/RR/temp/b1;	//粒界エネルギ密度（0.5J/m^2）を無次元化
	//delta=7.0;					//粒界幅（差分ブロック数にて表現）

	//K1=8.0*delta*gamma/PI/PI;	//勾配エネルギー係数[式(4.40)]
	K1=K0*delta*gamma/PI/PI;	//勾配エネルギー係数[式(4.40)]
	//W1=4.0*gamma/delta;				//ペナルティー項の係数[式(4.40)]
	W1=W0*gamma/delta;				//ペナルティー項の係数[式(4.40)]
	//amobi=1.0;
	//M1=amobi*PI*PI/(8.0*delta);	//粒界の易動度[式(4.40)]
	M1=amobi*PI*PI/(M0*delta);	//粒界の易動度[式(4.40)]
	//E1=50.0/RR/temp;				//粒界移動の駆動力
	E1=E0/RR/temp;				//粒界移動の駆動力

	time1=0.;					//計算カウント数の初期値
	time1max=10000001.0;		//計算カウント数の最大値

//*** 式(4.32) - 式(4.35)の配列（K,W,M,E）の設定 *****************************
	for(i=1;i<=nm;i++){
		for(j=1;j<=nm;j++){
			//wij[i][j]=W1;  aij[i][j]=K1;  tij[i][j]=0.0;  eij[i][j]=0.0;
			//if( (i==nm)||(j==nm) ){eij[i][j]=E1; tij[i][j]=M1;}
			//if(i>j){eij[i][j]=-eij[i][j];}
			//if(i==j){	wij[i][j]=0.0;  aij[i][j]=0.0;  tij[i][j]=0.0; eij[i][j]=0.0;}
			wij[i*ND+j]=W1;  aij[i*ND+j]=K1;  tij[i*ND+j]=0.0;  eij[i*ND+j]=0.0;
			if( (i==nm)||(j==nm) ){eij[i*ND+j]=E1; tij[i*ND+j]=M1;}
			if(i>j){eij[i*ND+j]=-eij[i*ND+j];}
			if(i==j){	wij[i*ND+j]=0.0;  aij[i*ND+j]=0.0;  tij[i*ND+j]=0.0; eij[i*ND+j]=0.0;}
		}
	}

//*** 初期場の設定と描画Window表示 *****************************************
	for(k=1;k<=nm;k++){
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				//ph[k][i][j]=0.0;  ph2[k][i][j]=0.0;//フェーズフィールド、および補助配列の初期化
				 ph[k*ND*ND+i*ND+j]=0.0;	//フェーズフィールドの初期化
				ph2[k*ND*ND+i*ND+j]=0.0;	//補助配列の初期化
			}
		}
	}
	ini000(ph, ND, N);	//初期場の設定

//**** シミュレーションスタート ******************************
//Nstep = 20;
iout = -1;
start: ;

	if((((int)(time1) % Nstep)==0)) {datsave(ph, ND, N);}	//一定繰返しカウント毎に場を保存
	if((((int)(time1) % Nstep)==0)) {datsave_paraview(ph, ND, N);}	//一定繰返しカウント毎に場を保存
	//if(time1==200.) {datsave();}					//特定の時間の場を保存

//**** 各差分プロックにおけるn00h[i][j]とm00h[n00][i][j]を調査 *********************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ip=i+1; im=i-1;  jp=j+1;  jm=j-1;
			if(i==ndm){ip=0;}  if(i==0){im=ndm;}//周期的境界条件
			if(j==ndm){jp=0;}  if(j==0){jm=ndm;}

//--- 位置(i,j)およびその周囲(i±1,j±1)において、pが０ではない方位の個数---
			n00=0;
			for(ii=1;ii<=nm;ii++){
				//if( (ph[ii][i][j]>0.0)||
				//  ( (ph[ii][i][j]==0.0)&&(ph[ii][ip][j]>0.0)||
				//  	(ph[ii][im][j]>0.0)||
				//  	(ph[ii][i][jp]>0.0)||
				//  	(ph[ii][i][jm]>0.0)   ) ){
				//  		n00++; m00h[n00][i][j]=ii;
				//}
				if( (ph[ii*ND*ND+i*ND+j]>0.0)||
					( (ph[ii*ND*ND+i*ND+j]==0.0)&&(ph[ii*ND*ND+ip*ND+j]>0.0)||
					  (ph[ii*ND*ND+im*ND+j]>0.0)||
					  (ph[ii*ND*ND+i*ND+jp]>0.0)||
					  (ph[ii*ND*ND+i*ND+jm]>0.0) ) ){
					  	n00++; m00h[n00*ND*ND+i*ND+j]=ii;
				}
			}
			//n00h[i][j]=n00;
			n00h[i*ND+j]=n00;
//--------------------------------------------------------------------------
		}
	}

//***** 発展方程式の計算 **********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ip=i+1; im=i-1;  jp=j+1;  jm=j-1;
			if(i==ndm){ip=0;}  if(i==0){im=ndm;}//周期的境界条件
			if(j==ndm){jp=0;}  if(j==0){jm=ndm;}

			//for(n1=1; n1<=n00h[i][j]; n1++){
			for(n1=1; n1<=n00h[i*ND+j]; n1++){
				//ii=m00h[n1][i][j];  pddtt=0.0;
				ii=m00h[n1*ND*ND+i*ND+j];  pddtt=0.0;
				//for(n2=1; n2<=n00h[i][j]; n2++){
				for(n2=1; n2<=n00h[i*ND+j]; n2++){
					//jj=m00h[n2][i][j];  sum1=0.0;
					jj=m00h[n2*ND*ND+i*ND+j];  sum1=0.0;
					//for(n3=1; n3<=n00h[i][j]; n3++){
					for(n3=1; n3<=n00h[i*ND+j]; n3++){
						//kk=m00h[n3][i][j];
						kk=m00h[n3*ND*ND+i*ND+j];
						//sum1+=0.5*(aij[ii][kk]-aij[jj][kk])*(ph[kk][ip][j]+ph[kk][im][j]
						//									+ph[kk][i][jp]+ph[kk][i][jm]-4.0*ph[kk][i][j])
						//								  +(wij[ii][kk]-wij[jj][kk])*ph[kk][i][j];	//[式(4.31)の一部]
						sum1+=0.5*(aij[ii*ND+kk]-aij[jj*ND+kk])*(ph[kk*ND*ND+ip*ND+j]+ph[kk*ND*ND+im*ND+j]
																+ph[kk*ND*ND+i*ND+jp]+ph[kk*ND*ND+i*ND+jm]-4.0*ph[kk*ND*ND+i*ND+j])
								 +(wij[ii*ND+kk]-wij[jj*ND+kk])*(ph[kk*ND*ND+i*ND+j]);	//[式(4.31)の一部]
					}
					//pddtt+=-2.0*tij[ii][jj]/double(n00h[i][j])
					//	*(sum1-8.0/PI*eij[ii][jj]*sqrt(ph[ii][i][j]*ph[jj][i][j]));
					pddtt+=-2.0*tij[ii*ND+jj]/double(n00h[i*ND+j])
				  *(sum1-8.0/PI*eij[ii*ND+jj]*sqrt(ph[ii*ND*ND+i*ND+j]*ph[jj*ND*ND+i*ND+j]));
					//フェーズフィールドの発展方程式[式(4.31)]
				}
				//ph2[ii][i][j]=ph[ii][i][j]+pddtt*delt;		//フェーズフィールドの時間発展（陽解法）
				//if(ph2[ii][i][j]>=1.0){ph2[ii][i][j]=1.0;}//フェーズフィールドの変域補正
				//if(ph2[ii][i][j]<=0.0){ph2[ii][i][j]=0.0;}
				ph2[ii*ND*ND+i*ND+j]=ph[ii*ND*ND+i*ND+j]+pddtt*delt;		//フェーズフィールドの時間発展（陽解法）
				if(ph2[ii*ND*ND+i*ND+j]>=1.0){ph2[ii*ND*ND+i*ND+j]=1.0;}//フェーズフィールドの変域補正
				if(ph2[ii*ND*ND+i*ND+j]<=0.0){ph2[ii*ND*ND+i*ND+j]=0.0;}
			}
		}//j
	}//i

	for(k=1;k<=nm;k++){
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				//ph[k][i][j]=ph2[k][i][j];	//補助配列を主配列に移動
				ph[k*ND*ND+i*ND+j]=ph2[k*ND*ND+i*ND+j];	//補助配列を主配列に移動
			}
		}
	}

//*** フェーズフィールドの規格化補正 ***********************************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//sum1=0.0; for(k=1;k<=nm;k++){ sum1+=ph[k][i][j]; }
			//for(k=1;k<=nm;k++){ ph[k][i][j]=ph[k][i][j]/sum1; }
			sum1=0.0; for(k=1;k<=nm;k++){ sum1+=ph[k*ND*ND+i*ND+j]; }
			for(k=1;k<=nm;k++){ ph[k*ND*ND+i*ND+j]=ph[k*ND*ND+i*ND+j]/sum1; }
		}
	}

//*********************************************************************
	//if(keypress()){return 0;}//キー待ち状態
	time1=time1+1.;  if(time1<time1max){goto start;}//最大カウント数に到達したかどうかの判断

	end:;
  return 0;
}

//************ 初期場(フェーズフィールド)の設定サブルーチン *************
void ini000(double *ph, int ND, int N)
{
	int i, j, k, l, it;		//整数
	int ii, jj, kk;				//整数
	int ip, im, jp, jm;		//整数
	int x1, y1, x1h[10], y1h[10];//初期核の座標
	double sum1, t, r0, phi, r;
	//srand(time(NULL)); // 乱数初期化
	int nd=ND, ndm=ND-1, nm=N-1, nmm=N-2;

	x1h[1]=0.2*nd;   y1h[1]=0.2*nd;		//初期核１の座標設定
	x1h[2]=0.75*nd;  y1h[2]=0.4*nd;		//初期核２の座標設定
	x1h[3]=0.5*nd;   y1h[3]=0.75*nd;	//初期核３の座標設定

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//for(ii=1;ii<=nm-1;ii++){ph[ii][i][j]=0.0;}
			//ph[nm][i][j]=1.0;//nm番目のフェーズフィールドを１に初期化
			for(ii=1;ii<=nm-1;ii++){ph[ii*ND*ND+i*ND+j]=0.0;}
			ph[nm*ND*ND+i*ND+j]=1.0;//nm番目のフェーズフィールドを１に初期化
		}
	}

	r0=10.0;
	for(ii=1;ii<=nm-1;ii++){
		x1=x1h[ii]; y1=y1h[ii];
		//x1=nd*DRND(1); y1=nd*DRND(1);
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				r=sqrt( (double(i-x1))*(double(i-x1))+(double(j-y1))*(double(j-y1)) );
				//if(r<=r0){ ph[ii][i][j]=1.0;  ph[nm][i][j]=0.0; } //初期核位置のフェーズフィールドを設定
				if(r<=r0){
					ph[ii*ND*ND+i*ND+j]=1.0;
					ph[nm*ND*ND+i*ND+j]=0.0;
				} //初期核位置のフェーズフィールドを設定
			}
		}
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//sum1=0.0; for(k=1;k<=nm;k++){ sum1+=ph[k][i][j]; }
			//for(k=1;k<=nm;k++){ ph[k][i][j]=ph[k][i][j]/sum1; }//フェーズフィールドの規格化補正
			sum1=0.0; for(k=1;k<=nm;k++){ sum1+=ph[k*ND*ND+i*ND+j]; }
			for(k=1;k<=nm;k++){ ph[k*ND*ND+i*ND+j]=ph[k*ND*ND+i*ND+j]/sum1; }//フェーズフィールドの規格化補正
		}
	}

}

//************ データ保存サブルーチン *******************************
void datsave(double *ph, int ND, int N)
{
	FILE		*stream;		//ストリームのポインタ設定
	int 		i, j, k;		//整数
	int nd=ND, ndm=ND-1, nm=N-1, nmm=N-2;

	stream = fopen("test.dat", "a");	//書き込む先のファイルを追記方式でオープン
	fprintf(stream, "%e  \n", time1);	//計算カウント数の保存
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=nm;k++){
				//fprintf(stream, "%e   ", ph[k][i][j]);//フェーズフィールドの保存
				fprintf(stream, "%e   ", ph[k*ND*ND+i*ND+j]);//フェーズフィールドの保存
			}
		}
	}
	fprintf(stream, "\n");	//改行の書き込み
	fclose(stream);					//ファイルをクローズ
}

void datsave_paraview(double *ph, int ND, int N)
{
	FILE	*fp;
	char	fName[256];
	int 	i, j, k;
	int nd=ND, ndm=ND-1, nm=N-1, nmm=N-2;
	double *pht  = (double *)malloc(sizeof(double)*( ND*ND + ND ));
	
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			pht[i*ND+j]=0.0;
		}
	}
	
	iout = iout + 1;
	printf("paraview output no.%06d \n",iout);
	for(k=1;k<=nm;k++){
		sprintf(fName,"mpf0_N%03d_result%06d.vtk",k,iout);
		fp = fopen(fName, "w");
		fprintf(fp,"# vtk DataFile Version 3.0 \n");
		fprintf(fp,"output.vtk \n");
		fprintf(fp,"ASCII \n");
		fprintf(fp,"DATASET STRUCTURED_POINTS \n");
		fprintf(fp,"DIMENSIONS %5d %5d %5d \n",(ndm+1),(ndm+1),1);
		fprintf(fp,"ORIGIN 0.0 0.0 0.0 \n");
		fprintf(fp,"ASPECT_RATIO 1 1 1 \n");
		fprintf(fp,"POINT_DATA %16d \n",((ndm+1)*(ndm+1)*1));
		fprintf(fp,"SCALARS phase_field float \n");
		fprintf(fp,"LOOKUP_TABLE default \n");
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				//fprintf(stream, "%e   ", ph[k][i][j]);//Save phase field
				//fprintf(fp,"%10.6f\n", ph[k][i][j]);
				fprintf(fp,"%10.6f\n", ph[k*ND*ND+i*ND+j]);
				pht[i*ND+j]+=ph[k*ND*ND+i*ND+j]*float(k);
			}
		}
		fclose(fp);
	}
	
	for(k=0;k<=0;k++){
		sprintf(fName,"mpf0_N%03d_result%06d.vtk",k,iout);
		fp = fopen(fName, "w");
		fprintf(fp,"# vtk DataFile Version 3.0 \n");
		fprintf(fp,"output.vtk \n");
		fprintf(fp,"ASCII \n");
		fprintf(fp,"DATASET STRUCTURED_POINTS \n");
		fprintf(fp,"DIMENSIONS %5d %5d %5d \n",(ndm+1),(ndm+1),1);
		fprintf(fp,"ORIGIN 0.0 0.0 0.0 \n");
		fprintf(fp,"ASPECT_RATIO 1 1 1 \n");
		fprintf(fp,"POINT_DATA %16d \n",((ndm+1)*(ndm+1)*1));
		fprintf(fp,"SCALARS phase_field float \n");
		fprintf(fp,"LOOKUP_TABLE default \n");
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				fprintf(fp,"%10.6f\n", pht[i*ND+j]);
			}
		}
		fclose(fp);
	}
	free(pht);
}
//*********** データ入力サブルーチン **************************
void datin(double *ph, int ND, int N)
{
	FILE		*datin0;//ストリームのポインタ設定
	int 		i, j, k;//整数
	int nd=ND, ndm=ND-1, nm=N-1, nmm=N-2;

	datin0 = fopen("test.dat", "r");//読み込み元のファイルをオープン
	fscanf(datin0, "%lf", &time1);	//計算カウント数の読み込み
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=nm;k++){
				//fscanf(datin0, "%lf  ", &ph[k][i][j]);//フェーズフィールドの読み込み
				fscanf(datin0, "%lf  ", &ph[k*ND*ND+i*ND+j]);//フェーズフィールドの読み込み
			}
		}
	}
	fclose(datin0);//ファイルをクローズ

}