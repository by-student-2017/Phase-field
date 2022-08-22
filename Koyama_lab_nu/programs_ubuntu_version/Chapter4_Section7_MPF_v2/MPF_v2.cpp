//データを読み込みNo.1-5
//粒番号は1-5
//界面エネルギー密度に方位差依存性を考慮（できるがしていない）
//濃度場も考慮（２次式）、β相の結晶粒は１個 No.5

#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

#define DRND(x) ((double)(x)/RAND_MAX*rand())//乱数の設定

//#define ND 100				//差分計算における計算領域一辺の分割数
//#define N 6					//考慮する結晶方位の数＋１

	//int nd=ND, ndm=ND-1;	//計算領域の一辺の差分分割数(差分ブロック数), ND-1を定義
	//int nm=N-1, nmm=N-2;	//考慮する結晶方位の数、N-2（考慮する結晶方位の数－１）を定義
	double PI=3.141592;		//pi
	double RR=8.3145;		//ガス定数
	double time1;			//計算カウント数
	//double ph[N][ND][ND], ph2[N][ND][ND];	//フェーズフィールド、フェーズフィールド補助配列
	double c0;				//平均組成
	//double ch[ND][ND], ch2[ND][ND];	//濃度場、濃度場の補助配列
	//double aij[N][N];		//勾配エネルギー係数
	//double wij[N][N];		//ペナルティー項の係数
	//double tij[N][N];		//粒界の易動度
	//double eij[N][N];		//粒界移動の駆動力
	//int m00h[N][ND][ND];	//位置(i,j)およびその周囲(i±1,j±1)において、pが０ではない方位の番号
	//int n00h[ND][ND];		//位置(i,j)およびその周囲(i±1,j±1)において、pが０ではない方位の個数
	int Nstep, iout;

	void ini000(double *ph, double *ch, int ND, int N);	//初期場の設定サブルーチン
	void datsave(double *ph, double *ch, int ND, int N);//データ保存サブルーチン
	void datsave_paraview(double *ph, double *ch, int ND, int N);	//データ保存サブルーチン
	void datin(double *ph, int ND, int N);	//データ入力サブルーチン

//******* メインプログラム ******************************************
int main(void)
{
	int ND, N;
	int nd, ndm, nm, nmm;
	
	int i, j, k, l, ii, jj, kk, ll, it;	//整数
	int ip, im, jp, jm;									//整数
	int n1, n2, n3;											//整数
	int nalph;
	int n00;//位置(i,j)およびその周囲(i±1,j±1)において、pが０ではない方位の個数
	//int n000;//位置(i,j)において、pが０ではない方位の個数（n00>=n000）
	//double cah[ND][ND], cbh[ND][ND];//局所平衡（平行接線則）の濃度場
	//double sah[ND][ND], sbh[ND][ND];//α相とβ相の存在確率（SαとSβ）
	double time1max;//計算カウント数の最大値（計算終了カウント）
	double delt, L, b1;//時間きざみ、計算領域１辺の長さ、差分ブロック一辺の長さ
	double M1, M0;		//粒界の易動度
	double W1, W0;		//ペナルティー項の係数
	double K1, K0;		//勾配エネルギー係数
	double E1, E0;		//粒界移動の駆動力
	double temp;		//温度
	double sum1, sum2, sum3;//各種の和の作業変数
	double pddtt;		//フェーズフィールドの時間変化率

	double gamma, gamma0;//粒界エネルギ密度
	double delta;	//粒界幅（差分ブロック数にて表現）
	double amobi;	//粒界の易動度
	double vm0;		//モル体積

	double A_alph, A_beta, c_alph0, c_beta0;//化学的自由エネルギー内のパラメータ
	double A_alph0, A_beta0;		//化学的自由エネルギー内のパラメータ
	double c, sa, sb, Da, Db, D;	//濃度場、Sα、Sβ、拡散係数（Dα、Dβ）、重み付拡散係数D
	double gii, gjj, cii, cjj;		//化学的自由エネルギー、局所平衡濃度
	double dgdc;					//化学的自由エネルギーの濃度微分
	double cddtt;					//濃度の変化率
	double dc0;						//濃度場補正用の作業変数
	double dev1_a, dev1_b, dev2_a, dev2_b;//拡散方程式内における微分計算の際の作業変数


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
	nalph   = data[15];
	A_alph0 = data[16];
	c_alph0 = data[17];
	A_beta0 = data[18];
	c_beta0 = data[19];
	Da      = data[20];
	Db      = data[21];
	D       = data[22];
	printf("---------------------------------\n");
	//
	nd=ND;
	ndm=ND-1;
	nm=N-1;
	nmm=N-2;
	//
	double *ph   = (double *)malloc(sizeof(double)*( N*ND*ND + ND*ND + ND ));	//フェーズフィールド
	double *ph2  = (double *)malloc(sizeof(double)*( N*ND*ND + ND*ND + ND ));	//フェーズフィールド補助配列
	double *ch   = (double *)malloc(sizeof(double)*( ND*ND + ND ));
	double *ch2  = (double *)malloc(sizeof(double)*( ND*ND + ND ));
	double *aij  = (double *)malloc(sizeof(double)*( ND*ND + ND ));	//勾配エネルギー係数
	double *wij  = (double *)malloc(sizeof(double)*( ND*ND + ND ));	//ペナルティー項の係数
	double *tij  = (double *)malloc(sizeof(double)*( ND*ND + ND ));	//粒界の易動度
	double *eij  = (double *)malloc(sizeof(double)*( ND*ND + ND ));	//粒界移動の駆動力
	double *m00h = (double *)malloc(sizeof(double)*( N*ND*ND + ND*ND + ND ));	//位置(i,j)およびその周囲(i±1,j±1)において、pが0ではない方位の番号
	double *n00h = (double *)malloc(sizeof(double)*( ND*ND + ND ));	//位置(i,j)およびその周囲(i±1,j±1)において、pが0ではない方位の個数
	//
	double *cah  = (double *)malloc(sizeof(double)*( ND*ND + ND ));	//局所平衡（平行接線則）の濃度場
	double *cbh  = (double *)malloc(sizeof(double)*( ND*ND + ND ));	//局所平衡（平行接線則）の濃度場
	double *sah  = (double *)malloc(sizeof(double)*( ND*ND + ND ));	//α相とβ相の存在確率（SαとSβ）
	double *sbh  = (double *)malloc(sizeof(double)*( ND*ND + ND ));	//α相とβ相の存在確率（SαとSβ）
	//
	//printf("delt(0.2)=  "); scanf(" %lf",&delt);//時間刻みの入力	//	delt=0.2;
	//printf("c0(0.4)=  "); scanf(" %lf",&c0);//平均組成の入力	//c0=0.2;

	//temp=1000.0; 							//温度(K)
	//L=2000.0;									//計算領域の一辺の長さ(nm)
	b1=L/(double)ND*1.0e-9; 	//差分プロック１辺の長さ(m)

	//vm0=7.0e-6;								//モル体積
	//gamma=0.5*vm0/RR/temp/b1;	//粒界エネルギ密度（0.5J/m^2）を無次元化
	gamma=gamma0*vm0/RR/temp/b1;	//粒界エネルギ密度（0.5J/m^2）を無次元化
	//delta=5.0;								//粒界幅（差分ブロック数にて表現）

	//K1=8.0*delta*gamma/PI/PI;	//勾配エネルギー係数[式(3.23)]
	K1=K0*delta*gamma/PI/PI;	//勾配エネルギー係数[式(3.23)]
	//W1=4.0*gamma/delta;				//ペナルティー項の係数[式(3.23)]
	W1=W0*gamma/delta;				//ペナルティー項の係数[式(3.23)]
	//amobi=1.0;
	//M1=amobi*PI*PI/(8.0*delta);	//粒界の易動度[式(3.23)]
	M1=amobi*PI*PI/(M0*delta);	//粒界の易動度[式(3.23)]
	//E1=50.0/RR/temp;					//粒界移動の駆動力
	E1=E0/RR/temp;					//粒界移動の駆動力

	time1=0.0;
	//time1max=10000001.;//計算カウント数の初期値と最大値

	//nalph=4;	//α相の結晶方位の数（考慮しているα相の結晶粒の数）
	//A_alph=1.0e+02/RR/temp;//化学的自由エネルギー内のパラメータ
	A_alph=A_alph0/RR/temp;//化学的自由エネルギー内のパラメータ
	//c_alph0=0.1;
	//A_beta=1.0e+02/RR/temp;
	A_beta=A_beta0/RR/temp;
	//c_beta0=0.9;
	//Da=Db=D=1.0;	//各相の拡散係数は全て等しいと仮定

//*** 式(4.36)-式(4.39)の配列（K,W,M,E）の設定 **************************
	for(i=1;i<=nm;i++){
		for(j=i+1;j<=nm;j++){
			//wij[i][j]=wij[j][i]=W1;
			//aij[i][j]=aij[j][i]=K1;
			//tij[i][j]=tij[j][i]=M1;
			//eij[i][j]=eij[j][i]=0.0;
			wij[i*ND+j]=wij[j*ND+i]=W1;
			aij[i*ND+j]=aij[j*ND+i]=K1;
			tij[i*ND+j]=tij[j*ND+i]=M1;
			eij[i*ND+j]=eij[j*ND+i]=0.0;
		}
		//wij[i][i]=0.0;  aij[i][i]=0.0;  tij[i][i]=0.0;  eij[i][i]=0.0;
		wij[i*ND+i]=0.0; aij[i*ND+i]=0.0; tij[i*ND+i]=0.0; eij[i*ND+i]=0.0;
	}

//*** 初期場の設定と描画Window表示 *****************************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=1;k<=nm;k++){
				//ph[k][i][j]=0.0;  ph2[k][i][j]=0.0;//フェーズフィールド、および補助配列の初期化
				 ph[k*ND*ND+i*ND+j]=0.0;	//フェーズフィールドの初期化
				ph2[k*ND*ND+i*ND+j]=0.0;	//補助配列の初期化
			}
			//ch[i][j]=0.0;  ch2[i][j]=0.0;//濃度場、および補助配列の初期化
			 ch[i*ND+j]=0.0;	//濃度場の初期化
			ch2[i*ND+j]=0.0;	//補助配列の初期化
		}
	}

	ini000(ph, ch, ND, N);//初期場の設定
	//datin(ph, ND, N);//初期場の入力

//**** シミュレーションスタート ******************************
//Nstep = 200;
iout = -1;
start: ;

	if((((int)(time1) % Nstep)==0)) {datsave(ph, ch, ND, N);}//一定繰返しカウント毎に場を保存
	if((((int)(time1) % Nstep)==0)) {datsave_paraview(ph, ch, ND, N);}//一定繰返しカウント毎に場を保存
	//if(time1==200.) {datsave();}//特定の時間の場を保存

//******  局所平衡組成（cαとcβ）の計算[式(4.47)]  ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//c=ch[i][j];	//濃度場
			c=ch[i*ND+j];	//濃度場
			//sum1=0.0; for(ii=1;ii<=nalph;ii++){ sum1+=ph[ii][i][j]; }
			sum1=0.0; for(ii=1;ii<=nalph;ii++){ sum1+=ph[ii*ND*ND+i*ND+j]; }
			//sa=sah[i][j]=sum1;	//Sαの計算[式(4.46)]
			sa=sah[i*ND+j]=sum1;	//Sαの計算[式(4.46)]
			//sum1=0.0; for(ii=nalph+1;ii<=nm;ii++){ sum1+=ph[ii][i][j]; }
			sum1=0.0; for(ii=nalph+1;ii<=nm;ii++){ sum1+=ph[ii*ND*ND+i*ND+j]; }
			//sb=sbh[i][j]=sum1;	//Sβの計算[式(4.46)]
			sb=sbh[i*ND+j]=sum1;	//Sβの計算[式(4.46)]

			//局所平衡組成の計算[式(4.47)]
			//cah[i][j]=(A_beta*c+(A_alph*c_alph0-A_beta*c_beta0)*sb)/(A_beta*sa+A_alph*sb);
			cah[i*ND+j]=(A_beta*c+(A_alph*c_alph0-A_beta*c_beta0)*sb)/(A_beta*sa+A_alph*sb);
			//cbh[i][j]=(A_alph*c+(A_beta*c_beta0-A_alph*c_alph0)*sa)/(A_beta*sa+A_alph*sb);
			cbh[i*ND+j]=(A_alph*c+(A_beta*c_beta0-A_alph*c_alph0)*sa)/(A_beta*sa+A_alph*sb);
			//if(cah[i][j]>=1.0){cah[i][j]=1.0;}  if(cah[i][j]<=0.0){cah[i][j]=0.0;}//濃度場の変域補正
			if(cah[i*ND+j]>=1.0){cah[i*ND+j]=1.0;}  if(cah[i*ND+j]<=0.0){cah[i*ND+j]=0.0;}//濃度場の変域補正
			//if(cbh[i][j]>=1.0){cbh[i][j]=1.0;}  if(cbh[i][j]<=0.0){cbh[i][j]=0.0;}
			if(cbh[i*ND+j]>=1.0){cbh[i*ND+j]=1.0;}  if(cbh[i*ND+j]<=0.0){cbh[i*ND+j]=0.0;}
		}
	}

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
				//		//printf("%d  ", n00);
				//}
				if( (ph[ii*ND*ND+i*ND+j]>0.0)||
				  ( (ph[ii*ND*ND+i*ND+j]==0.0)&&(ph[ii*ND*ND+ip*ND+j]>0.0)||
				  	(ph[ii*ND*ND+im*ND+j]>0.0)||
				  	(ph[ii*ND*ND+i*ND+jp]>0.0)||
				  	(ph[ii*ND*ND+i*ND+jm]>0.0)   ) ){
				  		n00++; m00h[n00*ND*ND+i*ND+j]=ii;
						//printf("%d  ", n00);
				}
			}
			//n00h[i][j]=n00;
			n00h[i*ND+j]=n00;
//--------------------------------------------------------------------------
		}
	}


//***** フェーズフィールドの発展方程式の計算 **********************************
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
						//		 +(wij[ii][kk]-wij[jj][kk])*ph[kk][i][j];//[式(4.41)の一部]
						sum1+=0.5*(aij[ii*ND+kk]-aij[jj*ND+kk])*(ph[kk*ND*ND+ip*ND+j]+ph[kk*ND*ND+im*ND+j]
																+ph[kk*ND*ND+i*ND+jp]+ph[kk*ND*ND+i*ND+jm]-4.0*ph[kk*ND*ND+i*ND+j])
								 +(wij[ii*ND+kk]-wij[jj*ND+kk])*(ph[kk*ND*ND+i*ND+j]);//[式(4.41)の一部]
					}

					//if(ii<=nalph){gii=A_alph*(cah[i][j]-c_alph0)*(cah[i][j]-c_alph0); cii=cah[i][j]; }
					if(ii<=nalph){gii=A_alph*(cah[i*ND+j]-c_alph0)*(cah[i*ND+j]-c_alph0); cii=cah[i*ND+j]; }
					//else{gii=A_beta*(cbh[i][j]-c_beta0)*(cbh[i][j]-c_beta0); cii=cbh[i][j]; }//式(4.43)
					else{gii=A_beta*(cbh[i*ND+j]-c_beta0)*(cbh[i*ND+j]-c_beta0); cii=cbh[i*ND+j]; }//式(4.43)
					//if(jj<=nalph){gjj=A_alph*(cah[i][j]-c_alph0)*(cah[i][j]-c_alph0); cjj=cah[i][j]; }
					if(jj<=nalph){gjj=A_alph*(cah[i*ND+j]-c_alph0)*(cah[i*ND+j]-c_alph0); cjj=cah[i*ND+j]; }
					//else{gjj=A_beta*(cbh[i][j]-c_beta0)*(cbh[i][j]-c_beta0); cjj=cbh[i][j]; }//式(4.43)
					else{gjj=A_beta*(cbh[i*ND+j]-c_beta0)*(cbh[i*ND+j]-c_beta0); cjj=cbh[i*ND+j]; }//式(4.43)
						//dgdc=2.0*A_alph*(cah[i][j]-c_alph0);//化学的自由エネルギーの濃度微分[式(4.48)]の計算
						dgdc=2.0*A_alph*(cah[i*ND+j]-c_alph0);//化学的自由エネルギーの濃度微分[式(4.48)]の計算
						//dgdc=2.0*A_beta*(cbh[i][j]-c_beta0);
						
						//pddtt+=-2.0*tij[ii][jj]/double(n00h[i][j])*(sum1+gii-gjj-(cii-cjj)*dgdc );
						pddtt+=-2.0*tij[ii*ND+jj]/double(n00h[i*ND+j])*(sum1+gii-gjj-(cii-cjj)*dgdc );
						//フェーズフィールドの発展方程式[式(4.41)]
				}
				//ph2[ii][i][j]=ph[ii][i][j]+pddtt*delt;//フェーズフィールドの時間発展（陽解法）
				//if(ph2[ii][i][j]>=1.0){ph2[ii][i][j]=1.0;}//フェーズフィールドの変域補正
				//if(ph2[ii][i][j]<=0.0){ph2[ii][i][j]=0.0;}
				ph2[ii*ND*ND+i*ND+j]=ph[ii*ND*ND+i*ND+j]+pddtt*delt;//フェーズフィールドの時間発展（陽解法）
				if(ph2[ii*ND*ND+i*ND+j]>=1.0){ph2[ii*ND*ND+i*ND+j]=1.0;}//フェーズフィールドの変域補正
				if(ph2[ii*ND*ND+i*ND+j]<=0.0){ph2[ii*ND*ND+i*ND+j]=0.0;}
			}
		}//j
	}//i

//*****  濃度場の時間発展の計算（拡散方程式）**********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ip=i+1;  im=i-1;  jp=j+1;  jm=j-1;
			if(i==ndm){ip=0;}  if(i==0){im=ndm;}//周期的境界条件
			if(j==ndm){jp=0;}  if(j==0){jm=ndm;}

			//拡散方程式内における微分計算
			//dev1_a=0.25*( (sah[ip][j]-sah[im][j])*(cah[ip][j]-cah[im][j])
			//			 +(sah[i][jp]-sah[i][jm])*(cah[i][jp]-cah[i][jm]) );
			dev1_a=0.25*( (sah[ip*ND+j]-sah[im*ND+j])*(cah[ip*ND+j]-cah[im*ND+j])
						 +(sah[i*ND+jp]-sah[i*ND+jm])*(cah[i*ND+jp]-cah[i*ND+jm]) );
			//dev1_b=0.25*( (sbh[ip][j]-sbh[im][j])*(cbh[ip][j]-cbh[im][j])
			//			 +(sbh[i][jp]-sbh[i][jm])*(cbh[i][jp]-cbh[i][jm]) );
			dev1_b=0.25*( (sbh[ip*ND+j]-sbh[im*ND+j])*(cbh[ip*ND+j]-cbh[im*ND+j])
						 +(sbh[i*ND+jp]-sbh[i*ND+jm])*(cbh[i*ND+jp]-cbh[i*ND+jm]) );
			//dev2_a=sah[i][j]*(cah[ip][j]+cah[im][j]+cah[i][jp]+cah[i][jm]-4.0*cah[i][j]);
			dev2_a=sah[i*ND+j]*(cah[ip*ND+j]+cah[im*ND+j]+cah[i*ND+jp]+cah[i*ND+jm]-4.0*cah[i*ND+j]);
			//dev2_b=sbh[i][j]*(cbh[ip][j]+cbh[im][j]+cbh[i][jp]+cbh[i][jm]-4.0*cbh[i][j]);
			dev2_b=sbh[i*ND+j]*(cbh[ip*ND+j]+cbh[im*ND+j]+cbh[i*ND+jp]+cbh[i*ND+jm]-4.0*cbh[i*ND+j]);

			cddtt=Da*(dev1_a+dev2_a)+Db*(dev1_b+dev2_b);	//拡散方程式[式(4.42)]
			//ch2[i][j]=ch[i][j]+cddtt*delt;	//濃度場の時間発展(陽解法)
			//ch2[i][j]=ch[i][j]+cddtt*delt+(2.*DRND(1.)-1.)*0.001;	//濃度場の時間発展(陽解法)
			ch2[i*ND+j]=ch[i*ND+j]+cddtt*delt+(2.*DRND(1.)-1.)*0.001;	//濃度場の時間発展(陽解法)
		}
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=1;k<=nm;k++){
				//ph[k][i][j]=ph2[k][i][j];//補助配列を主配列に移動（フェーズフィールド）
				ph[k*ND*ND+i*ND+j]=ph2[k*ND*ND+i*ND+j];//補助配列を主配列に移動（フェーズフィールド）
			}
			//ch[i][j]=ch2[i][j];//補助配列を主配列に移動（濃度場）
			ch[i*ND+j]=ch2[i*ND+j];//補助配列を主配列に移動（濃度場）
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

//*** 濃度場の収支補正 *************************************************************
	//sum1=0.0; for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ sum1+=ch[i][j]; } }
	sum1=0.0; for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ sum1+=ch[i*ND+j]; } }
	dc0=sum1/nd/nd-c0;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//ch[i][j]=ch[i][j]-dc0;
			//if(ch[i][j]>1.0){ch[i][j]=1.0;}   if(ch[i][j]<0.0){ch[i][j]=0.0;}
			ch[i*ND+j]=ch[i*ND+j]-dc0;
			if(ch[i*ND+j]>1.0){ch[i*ND+j]=1.0;}
			if(ch[i*ND+j]<0.0){ch[i*ND+j]=0.0;}
	  }
	}

//*********************************************************************
	//if(keypress()){return 0;}//キー待ち状態
	time1=time1+1.;  if(time1<time1max){goto start;}//最大カウント数に到達したかどうかの判断

	end:;
  return 0;
}

//************ 初期場(フェーズフィールドと濃度場)の設定サブルーチン *************
void ini000(double *ph, double *ch, int ND, int N)
{
	int i, j, k, l, it;		//整数
	int ii, jj, kk;				//整数
	int ip, im, jp, jm;		//整数
	int igx, igy, ixmin, ixmax, iymin, iymax;	//差分ブロック座標系の設定
	double x, y, xmin, xmax, ymin, ymax;			//規格化座標系の設定
	double sum1, sum2, t, r0, phi, r;					//作業変数
	//srand(time(NULL)); //乱数初期化
	int nd=ND, ndm=ND-1, nm=N-1, nmm=N-2;

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//ch[i][j]=0.5;//α相の濃度場（過飽和固溶体）を設定
			ch[i*ND+j]=0.5;//α相の濃度場（過飽和固溶体）を設定
		}
	}

	datin(ph, ND, N);//α相の結晶組織データを読み込み

	xmin=-1.0; xmax=1.0; ymin=-1.0;  ymax=1.0;	//規格化座標系
	ixmin=0; ixmax=ND; iymin=0;  iymax=ND;			//差分ブロック座標系

	r0=2.0;//β相のサイズ

//以下、α相の粒界の３重点（６箇所）にβ相の核を置く
	x=0.5;           igx=(ixmax-ixmin)/(xmax-xmin)*(x-xmin)+ixmin;
	y=0.5/sqrt(3.);  igy=(iymax-iymin)/(ymax-ymin)*(y-ymin)+iymin;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			r=sqrt( (double(i-igx))*(double(i-igx))+(double(j-igy))*(double(j-igy)) );
			//if(r<=r0){ ph[nm][i][j]=1.0; 	for(k=1;j<=nm-1;k++){ph[k][i][j]=0.0;} ch[i][j]=0.9; }
			if(r<=r0){ ph[nm*ND*ND+i*ND+j]=1.0;
				for(k=1;j<=nm-1;k++){ph[k*ND*ND+i*ND+j]=0.0;} 
				ch[i*ND+j]=0.9;
			}
		}
	}

	x=-0.5;          igx=(ixmax-ixmin)/(xmax-xmin)*(x-xmin)+ixmin;
	y=0.5/sqrt(3.);  igy=(iymax-iymin)/(ymax-ymin)*(y-ymin)+iymin;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			r=sqrt( (double(i-igx))*(double(i-igx))+(double(j-igy))*(double(j-igy)) );
			//if(r<=r0){ ph[nm][i][j]=1.0; 	for(k=1;j<=nm-1;k++){ph[k][i][j]=0.0;} ch[i][j]=0.9;  }
			if(r<=r0){
				ph[nm*ND*ND+i*ND+j]=1.0;
				for(k=1;j<=nm-1;k++){ph[k*ND*ND+i*ND+j]=0.0;}
				ch[i*ND+j]=0.9;
			}
		}
	}

	x=0.5;           igx=(ixmax-ixmin)/(xmax-xmin)*(x-xmin)+ixmin;
	y=-0.5/sqrt(3.); igy=(iymax-iymin)/(ymax-ymin)*(y-ymin)+iymin;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			r=sqrt( (double(i-igx))*(double(i-igx))+(double(j-igy))*(double(j-igy)) );
			//if(r<=r0){ ph[nm][i][j]=1.0; 	for(k=1;j<=nm-1;k++){ph[k][i][j]=0.0;} ch[i][j]=0.9;  }
			if(r<=r0){
				ph[nm*ND*ND+i*ND+j]=1.0;
				for(k=1;j<=nm-1;k++){ph[k*ND*ND+i*ND+j]=0.0;}
				ch[i*ND+j]=0.9;
			}
		}
	}

	x=-0.5;           igx=(ixmax-ixmin)/(xmax-xmin)*(x-xmin)+ixmin;
	y=-0.5/sqrt(3.);  igy=(iymax-iymin)/(ymax-ymin)*(y-ymin)+iymin;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			r=sqrt( (double(i-igx))*(double(i-igx))+(double(j-igy))*(double(j-igy)) );
			//if(r<=r0){ ph[nm][i][j]=1.0; 	for(k=1;j<=nm-1;k++){ph[k][i][j]=0.0;} ch[i][j]=0.9;  }
			if(r<=r0){
				ph[nm*ND*ND+i*ND+j]=1.0;
				for(k=1;j<=nm-1;k++){ph[k*ND*ND+i*ND+j]=0.0;}
				ch[i*ND+j]=0.9;
			}
		}
	}

	x=0.0;          igx=(ixmax-ixmin)/(xmax-xmin)*(x-xmin)+ixmin;
	y=1./sqrt(3.);  igy=(iymax-iymin)/(ymax-ymin)*(y-ymin)+iymin;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			r=sqrt( (double(i-igx))*(double(i-igx))+(double(j-igy))*(double(j-igy)) );
			//if(r<=r0){ ph[nm][i][j]=1.0; 	for(k=1;j<=nm-1;k++){ph[k][i][j]=0.0;} ch[i][j]=0.9;  }
			if(r<=r0){
				ph[nm*ND*ND+i*ND+j]=1.0;
				for(k=1;j<=nm-1;k++){ph[k*ND*ND+i*ND+j]=0.0;}
				ch[i*ND+j]=0.9;
			}
		}
	}

	x=0.0;           igx=(ixmax-ixmin)/(xmax-xmin)*(x-xmin)+ixmin;
	y=-1./sqrt(3.);  igy=(iymax-iymin)/(ymax-ymin)*(y-ymin)+iymin;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			r=sqrt( (double(i-igx))*(double(i-igx))+(double(j-igy))*(double(j-igy)) );
			//if(r<=r0){ ph[nm][i][j]=1.0; 	for(k=1;j<=nm-1;k++){ph[k][i][j]=0.0;} ch[i][j]=0.9;  }
			if(r<=r0){
				ph[nm*ND*ND+i*ND+j]=1.0;
				for(k=1;j<=nm-1;k++){ph[k*ND*ND+i*ND+j]=0.0;}
				ch[i*ND+j]=0.9;
			}
		}
	}

//--- フェーズフィールドの規格化補正と平均組成の算出 ------------------------
	sum2=0.0;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//sum1=0.0; for(k=1;k<=nm;k++){ sum1+=ph[k][i][j]; }
			//for(k=1;k<=nm;k++){ ph[k][i][j]=ph[k][i][j]/sum1; }
			//sum2+=ch[i][j];
			sum1=0.0; for(k=1;k<=nm;k++){ sum1+=ph[k*ND*ND+i*ND+j]; }
			for(k=1;k<=nm;k++){ ph[k*ND*ND+i*ND+j]=ph[k*ND*ND+i*ND+j]/sum1; }
			sum2+=ch[i*ND+j];
		}
	}
	c0=sum2/nd/nd;

}

//************ データ保存サブルーチン *******************************
void datsave(double *ph, double *ch, int ND, int N)
{
	FILE		*stream;//ストリームのポインタ設定
	int 		i, j, k;//整数
	double 	col;
	int nd=ND, ndm=ND-1, nm=N-1, nmm=N-2;

	stream = fopen("test.dat", "a");//書き込む先のファイルを追記方式でオープン
	fprintf(stream, "%e  \n", time1);//計算カウント数の保存
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=nm;k++){
				//fprintf(stream, "%e   ", ph[k][i][j]);//フェーズフィールドの保存
				fprintf(stream, "%e   ", ph[k*ND*ND+i*ND+j]);//フェーズフィールドの保存
			}
			//fprintf(stream, "%e   ", ch[i][j]);//濃度場の保存
			fprintf(stream, "%e   ", ch[i*ND+j]);//濃度場の保存
		}
	}
	fprintf(stream, "\n");//改行の書き込み
	fclose(stream);//ファイルをクローズ
}

void datsave_paraview(double *ph, double *ch, int ND, int N)
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
	for(k=0;k<=nm;k++){
		sprintf(fName,"mpf_N%03d_result%06d.vtk",k,iout);
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
				//fprintf(stream, "%e   ", ph[k][i][j]);//フェーズフィールドの保存
				//fprintf(fp,"%10.6f\n", ph[k][i][j]);
				fprintf(fp,"%10.6f\n", ph[k*ND*ND+i*ND+j]);
				pht[i*ND+j]+=ph[k*ND*ND+i*ND+j]*float(k+1.0);
			}
		}
		fclose(fp);
	}
	
	sprintf(fName,"mpf_N%03d_result%06d.vtk",(nm+1),iout);
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
	fprintf(fp,"SCALARS concentration float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			//fprintf(stream, "%e   ", ch[i][j]);//濃度場の保存
			//fprintf(fp,"%10.6f\n", ch[i][j]);
			fprintf(fp,"%10.6f\n", ch[i*ND+j]);
		}
	}
	fclose(fp);
	free(pht);
}
//*********** データ入力サブルーチン **************************
void datin(double *ph, int ND, int N)
{
	FILE		*datin0;//ストリームのポインタ設定
	int 		i, j, k;//整数
	double dami1;//作業変数
	int nd=ND, ndm=ND-1, nm=N-1, nmm=N-2;

	datin0 = fopen("test.dat", "r");//読み込み元のファイルをオープン
	fscanf(datin0, "%lf  ", &dami1);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=1;k<=nm;k++){
				//fscanf(datin0, "%lf  ", &ph[k][i][j]);//フェーズフィールドの読み込み
				fscanf(datin0, "%lf  ", &ph[k*ND*ND+i*ND+j]);//フェーズフィールドの読み込み
			}
		}
	}
	fclose(datin0);//ファイルをクローズ

}

