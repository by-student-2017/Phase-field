#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

#define DRND(x) ((double)(x)/RAND_MAX*rand())//乱数の関数設定

//#define ND 128						//差分計算における計算領域一辺の分割数(高速フーリエ変換を用いるため２のべき乗)
//#define IG 7						// 2^IG=ND

	//int nd=ND, ndm=ND-1; 			//計算領域の一辺の差分分割数(差分ブロック数)、ND-1を定義
	//int nd2=ND/2;				 	//ND/2を定義：高速フ－リエ変換で使用
	//int ig=IG;						//2^ig=ND
	double PI=3.14159;				//円周率
	double RR=8.3145;				//ガス定数
	double time1;					//計算カウント数(時間に比例)
	double ep0=8.8541878e-12; 		//真空の誘電率(F/m)
	double Peq;						//分極モーメントの平衡値

	//double s1h[ND][ND], s2h[ND][ND];//x方向の分極モーメント,y方向の分極モーメント

	double qs;						//フ－リエ変換(qs:-1)とフ－リエ逆変換(qs:1)の区別
	//double xi[ND][ND], xr[ND][ND], xif[ND], xrf[ND];//フ－リエ変換の実部・虚部配列
	//double s[ND],c[ND];				//sinとcosのテーブル
	//int ik[ND];						//ビット反転テーブル
	int iout;

	void ini000_s12(double *s1h, double *s2h, int ND);			//時間0における分極モーメントの初期プロファイル
	void table(double *s, double *c, int *ik, int ND, int ig);	//sinとcosのテーブルとビット反転テーブルの作成サブル－チン
	void fft(double *xrf, double *xif, double *s, double *c, int ND, int ig);	//１次元高速フーリエ変換
	void rcfft(double *xrf, double *xif, double *xr, double *xi, int *ik, double *s, double *c, int ND, int ig);//２次元高速フーリエ変換
	void datsave(double *s1h, double *s2h, int ND);	//デ－タ保存サブル－チン
	void datsave_paraview(double *s1h, double *s2h, int ND);//デ－タ保存サブル－チン
	void datin(double *s1h, double *s2h, int ND);	//初期データ読み込み

//******* メインプログラム ******************************************
int main(void)
{
    int ND;
	int nd, ndm, nd2, ig;
	
	int   i, j, k, l, ii, jj, Nstep;		//整数
	int   ip, im, jp, jm;					//整数

	//double s1qrh[ND][ND], s1qih[ND][ND];	//s1のフーリエ変換（実部、虚部）
	//double s2qrh[ND][ND], s2qih[ND][ND];	//s2のフーリエ変換（実部、虚部）
	//double s1h2[ND][ND], s2h2[ND][ND];		//s1とs2の補助配列

	//double ss1qrh[ND][ND], ss1qih[ND][ND];	//s1*s1のフーリエ変換（実部、虚部）
	//double ss2qrh[ND][ND], ss2qih[ND][ND];	//s2*s2のフーリエ変換（実部、虚部）
	//double s1s2qrh[ND][ND], s1s2qih[ND][ND];//s1*s2のフーリエ変換（実部、虚部）
	double ss1ss2;							//変域に関する数値計算誤差の補正用の作業変数

	double a0_aa, a0_a, a0_c;				//BaTiO3の格子定数(正方晶)
	double al, temp, delt, vm0;				//計算領域１編の長さ、温度、時間きざみ、モル体積
	double time1max;						//計算カウント数の最大値（計算終了カウント）
	double b1;								//差分ブロック１辺の長さ

	double s1, s2;							//分極モーメントのxおよびy成分
	double s1ip, s1im, s1jp, s1jm;			//s1の左右上下の値
	double s2ip, s2im, s2jp, s2jm;			//s2の左右上下の値

	double s1k, s1k_chem, s1k_surf, s1k_str, s1k_ddi;	//s1に関するポテンシャル
	//double s1k_dd[ND][ND];					//双極子-双極子相互作用ポテンシャル
	double s2k, s2k_chem, s2k_surf, s2k_str, s2k_ddi;	//s2に関するポテンシャル
	//double s2k_dd[ND][ND];					//双極子-双極子相互作用ポテンシャル
	double smob1, smob2;					//ドメイン界面の易動度
	double s1ddtt, s2ddtt;					//発展方程式の左辺

	double A1, A11, A12;					//化学的自由エネルギー内のパラメータ
	double A1e, A11e, A12e;
	double A111, A112, A123;
	double A111e, A112e, A123e;
	double A1111, A1112, A1122, A1123;
	double A1111e, A1112e, A1122e, A1123e;
	double Tc0;

	double nx, ny, alnn;					//逆空間の単位ベクトル成分と逆格子ベクトル長
	double kapaP, kapaPc;					//勾配エネルギー係数

	double E1_ex, E1_ex_0; 					//外部電場
	double ep00;							//比誘電率
	double Add0, Add0c;						//双極子-双極子相互作用計算における係数

	double t1, t2, t3, tQ, tR, tS, tT;		//平衡分極モーメント算出時の作業変数

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
	delt    = data[1];
	time1max= int(data[2]);
	Nstep   = int(data[3]);
	temp    = data[4];
	al      = data[5];
	vm0     = data[6];
	smob1   = data[7];
	smob2   = data[8];
	A1e     = data[9];
	A11e    = data[10];
	A12e    = data[11];
	A111e   = data[12];
	A112e   = data[13];
	A123e   = data[14];
	A1111e  = data[15];
	A1112e  = data[16];
	A1122e  = data[17];
	A1123e  = data[18];
	Tc0     = data[19];
	kapaPc  = data[20];
	E1_ex_0 = data[21];
	ep00    = data[22];
	Add0c   = data[23];
	printf("---------------------------------\n");
	//
	ig = int(log2(ND));
	nd=ND;
	ndm=ND-1;
	nd2=ND/2;
	//
	double *s1h     = (double *)malloc(sizeof(double)*( ND*ND ));	//x方向の分極モーメント
	double *s2h     = (double *)malloc(sizeof(double)*( ND*ND ));	//y方向の分極モーメント
	//
	double *xi      = (double *)malloc(sizeof(double)*( ND*ND ));//フ－リエ変換の実部・虚部配列
	double *xr      = (double *)malloc(sizeof(double)*( ND*ND ));//フ－リエ変換の実部・虚部配列
	//
	double *xif     = (double *)malloc(sizeof(double)*( ND ));//フ－リエ変換の実部・虚部配列
	double *xrf     = (double *)malloc(sizeof(double)*( ND ));//フ－リエ変換の実部・虚部配列
	//
	double *s       = (double *)malloc(sizeof(double)*( ND ));//sinのテーブル
	double *c       = (double *)malloc(sizeof(double)*( ND ));//cosのテーブル
	//
	int *ik         = (int *)malloc(sizeof(int)*( ND ));//ビット反転テーブル
	//
	double *s1qrh   = (double *)malloc(sizeof(double)*( ND*ND ));	//s1のフーリエ変換（実部）
	double *s1qih   = (double *)malloc(sizeof(double)*( ND*ND ));	//s1のフーリエ変換（虚部）
	double *s2qrh   = (double *)malloc(sizeof(double)*( ND*ND ));	//s2のフーリエ変換（実部）
	double *s2qih   = (double *)malloc(sizeof(double)*( ND*ND ));	//s2のフーリエ変換（虚部）
	double *s1h2    = (double *)malloc(sizeof(double)*( ND*ND ));	//s1の補助配列
	double *s2h2    = (double *)malloc(sizeof(double)*( ND*ND ));	//s2の補助配列
	double *ss1qrh  = (double *)malloc(sizeof(double)*( ND*ND ));	//s1*s1のフーリエ変換（実部）
	double *ss1qih  = (double *)malloc(sizeof(double)*( ND*ND ));	//s1*s1のフーリエ変換（虚部）
	double *ss2qrh  = (double *)malloc(sizeof(double)*( ND*ND ));	//s2*s2のフーリエ変換（実部）
	double *ss2qih  = (double *)malloc(sizeof(double)*( ND*ND ));	//s2*s2のフーリエ変換（虚部）
	double *s1s2qrh = (double *)malloc(sizeof(double)*( ND*ND ));	//s1*s2のフーリエ変換（実部）
	double *s1s2qih = (double *)malloc(sizeof(double)*( ND*ND ));	//s1*s2のフーリエ変換（虚部）
	//
	double *s1k_dd  = (double *)malloc(sizeof(double)*( ND*ND ));	//双極子-双極子相互作用ポテンシャル
	double *s2k_dd  = (double *)malloc(sizeof(double)*( ND*ND ));	//双極子-双極子相互作用ポテンシャル
	//
	//printf("DELT(0.1)=  "); scanf(" %lf",&delt);//時間刻み設定	//delt=0.1;

	time1=0.0;								//初期計算カウント数の設定
	//time1max=1.0+1.0e+07;					//最大計算カウント数の設定
	//時間は全て無次元化されているので注意。

	//temp=298.0;								//温度(K)

	//al=250.;								//計算領域の１辺の長さ（μm）
	b1=al*1.0E-06/nd;						//差分ブロック１辺の長さ（m）

	//a0_aa=0.4;  a0_a=0.3992;  a0_c=0.4036;	//格子定数(nm)
	//vm0=a0_a*a0_a*a0_c*1.0e-27*6.02e+23;	//モル体積（分子１モル）

	//smob1=1.; smob2=1.;						//構造相転移における易動度（規格化して１に設定）

//--- 化学的自由エネルギー内のパラメータ値 [表4.7参照]----------------------
	//A1=4.124e+05*vm0/RR/temp;
	A1=A1e*vm0/RR/temp;
	//A11=-2.097e+08*vm0/RR/temp;
	A11=A11e*vm0/RR/temp;
	//A12=7.974e+08*vm0/RR/temp;
	A12=A12e*vm0/RR/temp;
	//A111=1.294e+09*vm0/RR/temp;
	A111=A111e*vm0/RR/temp;
	//A112=-1.950e+09*vm0/RR/temp;
	A112=A112e*vm0/RR/temp;
	//A123=-2.500e+09*vm0/RR/temp;
	A123=A123e*vm0/RR/temp;
	//A1111=3.863e+10*vm0/RR/temp;
	A1111=A1111e*vm0/RR/temp;
	//A1112=2.529e+10*vm0/RR/temp;
	A1112=A1112e*vm0/RR/temp;
	//A1122=1.637e+10*vm0/RR/temp;
	A1122=A1122e*vm0/RR/temp;
	//A1123=1.367e+10*vm0/RR/temp;
	A1123=A1123e*vm0/RR/temp;

	//Tc0=115.0+273.0;  //K

//--- 平衡分極モーメント値の計算 ----------------------
	t1=3.0*A111/4.0/A1111;
	t2=A11/2.0/A1111;
	t3=A1*(temp-Tc0)/4.0/A1111;
	tQ=(3.0*t2-t1*t1)/9.0;
	tR=(9.0*t1*t2-27.0*t3-2.0*t1*t1*t1)/54.0;
	tS=pow((tR+sqrt(tQ*tQ*tQ+tR*tR)),(1.0/3.0));
	tT=pow((tR-sqrt(tQ*tQ*tQ+tR*tR)),(1.0/3.0));
	Peq=sqrt(fabs(tS+tR-t1/3.));//平衡分極モーメント値
	printf("Peq=  %f  \n", Peq);//平衡分極モーメント値の表示

	//kapaP=1.0e-15/ep0*vm0/RR/temp/b1/b1;//勾配エネルギー係数
	kapaP=kapaPc/ep0*vm0/RR/temp/b1/b1;//勾配エネルギー係数

	//E1_ex_0=0.0;//外部電場(x方向)
	//E1_ex_0=2.0e+6*vm0/RR/temp;//外部電場(x方向)

	//ep00=200.0;//比誘電率

	//Add0=1.0/ep0/ep00*vm0/RR/temp;//双極子-双極子相互作用計算における係数[式(4.52)参照]
	Add0=Add0c/ep0/ep00*vm0/RR/temp;//双極子-双極子相互作用計算における係数

//*** sinおよびcosテ－ブル、ビット反転テーブル、および初期場の設定 ***************

	table(s, c, ik, ND, ig);		//sinおよびcosテ－ブルとビット反転テーブルの設定
	ini000_s12(s1h, s2h, ND);		//時間0sにおける分極モーメントの初期プロファイル
	//datin(s1h, s2h, ND);			//初期組織場の入力

//**** シミュレーションスタート ******************************
//Nstep = 10;
start: ;

	//if(time1<=100.){Nstep=10;} else{Nstep=100;}//データ保存する時間間隔の変更
	if((((int)(time1) % Nstep)==0)){datsave(s1h, s2h, ND);} //一定繰返しカウント毎に組織データを保存
	if((((int)(time1) % Nstep)==0)){datsave_paraview(s1h, s2h, ND);} //一定繰返しカウント毎に組織データを保存
	//if((int)(time1)==2000){datsave(s1h, s2h, ND);} //特定計算カウントにおけるデ－タの保存

//**** s1のフーリエ変換[式(3.37)] ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//xr[i][j]=s1h[i][j]; xi[i][j]=0.;
			xr[i*ND+j]=s1h[i*ND+j];
			xi[i*ND+j]=0.;
		}
	}
	qs=-1.; rcfft(xrf, xif, xr, xi, ik, s, c, ND, ig);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//s1qrh[i][j]=xr[i][j];  s1qih[i][j]=xi[i][j];
			s1qrh[i*ND+j]=xr[i*ND+j];
			s1qih[i*ND+j]=xi[i*ND+j];
		}
	}
	//s1qrh[0][0]=s1qih[0][0]=0.;
	s1qrh[0]=s1qih[0]=0.0;

//**** s2のフーリエ変換[式(3.37)] ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//xr[i][j]=s2h[i][j]; xi[i][j]=0.;
			xr[i*ND+j]=s2h[i*ND+j];
			xi[i*ND+j]=0.0;
		}
	}
	qs=-1.; rcfft(xrf, xif, xr, xi, ik, s, c, ND, ig);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//s2qrh[i][j]=xr[i][j];  s2qih[i][j]=xi[i][j];
			s2qrh[i*ND+j]=xr[i*ND+j];
			s2qih[i*ND+j]=xi[i*ND+j];
		}
	}
	//s2qrh[0][0]=s2qih[0][0]=0.;
	s2qrh[0]=s2qih[0]=0.0;

//***** s1の双極子-双極子相互作用の計算[式(3.33)内の畳込み積分] ***********************
	for(i=0;i<=ndm;i++){
		if(i<=nd2-1){ii=i;}  if(i>=nd2){ii=i-nd;}
		for(j=0;j<=ndm;j++){
			if(j<=nd2-1){jj=j;}  if(j>=nd2){jj=j-nd;}
			alnn=sqrt((double)ii*(double)ii+(double)jj*(double)jj);
			if(alnn==0.){alnn=1.;}
			nx=(double)ii/alnn;  ny=(double)jj/alnn;
			//xr[i][j]=s1qrh[i][j]*nx*nx+s2qrh[i][j]*nx*ny;
			//xi[i][j]=s1qih[i][j]*nx*nx+s2qih[i][j]*nx*ny;
			xr[i*ND+j]=s1qrh[i*ND+j]*nx*nx+s2qrh[i*ND+j]*nx*ny;
			xi[i*ND+j]=s1qih[i*ND+j]*nx*nx+s2qih[i*ND+j]*nx*ny;
		}
	}
	qs=1.; rcfft(xrf, xif, xr, xi, ik, s, c, ND, ig);
	//for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ s1k_dd[i][j]=xr[i][j]; } }
	for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ s1k_dd[i*ND+j]=xr[i*ND+j]; } }

//***** s2の双極子-双極子相互作用の計算[式(3.33)内の畳込み積分] ************************
	for(i=0;i<=ndm;i++){
		if(i<=nd2-1){ii=i;}  if(i>=nd2){ii=i-nd;}
		for(j=0;j<=ndm;j++){
			if(j<=nd2-1){jj=j;}  if(j>=nd2){jj=j-nd;}
			alnn=sqrt((double)ii*(double)ii+(double)jj*(double)jj);
			if(alnn==0.){alnn=1.;}
			nx=(double)ii/alnn;  ny=(double)jj/alnn;
			//xr[i][j]=s1qrh[i][j]*ny*nx+s2qrh[i][j]*ny*ny;
			//xi[i][j]=s1qih[i][j]*ny*nx+s2qih[i][j]*ny*ny;
			xr[i*ND+j]=s1qrh[i*ND+j]*ny*nx+s2qrh[i*ND+j]*ny*ny;
			xi[i*ND+j]=s1qih[i*ND+j]*ny*nx+s2qih[i*ND+j]*ny*ny;
		}
	}
	qs=1.; rcfft(xrf, xif, xr, xi, ik, s, c, ND, ig);
	//for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ s2k_dd[i][j]=xr[i][j]; } }
	for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ s2k_dd[i*ND+j]=xr[i*ND+j]; } }

//******  非線形発展方程式の数値計算  ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ip=i+1;  im=i-1;  jp=j+1;  jm=j-1;
			if(i==ndm){ip=0;}  if(i==0){im=ndm;}//周期的境界条件
			if(j==ndm){jp=0;}  if(j==0){jm=ndm;}
			//s1=s1h[i][j]; s1ip=s1h[ip][j]; s1im=s1h[im][j]; s1jp=s1h[i][jp]; s1jm=s1h[i][jm];
			//s2=s2h[i][j]; s2ip=s2h[ip][j]; s2im=s2h[im][j]; s2jp=s2h[i][jp]; s2jm=s2h[i][jm];
			s1=s1h[i*ND+j]; s1ip=s1h[ip*ND+j]; s1im=s1h[im*ND+j]; s1jp=s1h[i*ND+jp]; s1jm=s1h[i*ND+jm];
			s2=s2h[i*ND+j]; s2ip=s2h[ip*ND+j]; s2im=s2h[im*ND+j]; s2jp=s2h[i*ND+jp]; s2jm=s2h[i*ND+jm];

			//勾配ポテンシャルの計算
			s1k_surf=-kapaP*(s1ip+s1im+s1jp+s1jm-4.*s1);
			s2k_surf=-kapaP*(s2ip+s2im+s2jp+s2jm-4.*s2);

			//化学ポテンシャルの計算[式(4.55)]
			s1k_chem=2.0*A1*(temp-Tc0)*s1
							+4.0*A11*s1*s1*s1+2.0*A12*s1*s2*s2
							+6.0*A111*pow(s1,5.0)+A112*(2.0*s1*pow(s2,4.0)+4.0*s2*s2*pow(s1,3.0))
							+8.0*A1111*pow(s1,7.0)
							+A1112*(6.0*pow(s1,5.0)*s2*s2+2.0*pow(s2,6.0)*s1)
							+4.0*A1122*pow(s1,3.0)*pow(s2,4.0);

			s2k_chem=2.0*A1*(temp-Tc0)*s2
							+4.0*A11*s2*s2*s2+2.0*A12*s2*s1*s1
							+6.0*A111*pow(s2,5.0)+A112*(2.0*s2*pow(s1,4.0)+4.0*s1*s1*pow(s2,3.0))
							+8.0*A1111*pow(s2,7.0)
							+A1112*(6.0*pow(s2,5.0)*s1*s1+2.0*pow(s1,6.0)*s2)
							+4.0*A1122*pow(s2,3.0)*pow(s1,4.0);

			//双極子-双極子ポテンシャルの計算[式(4.56)右辺第１項] 
			//s1k_ddi=Add0*s1k_dd[i][j];
			//s2k_ddi=Add0*s2k_dd[i][j];
			s1k_ddi=Add0*s1k_dd[i*ND+j];
			s2k_ddi=Add0*s2k_dd[i*ND+j];

			//全ポテンシャルの計算
			s1k=s1k_chem+s1k_surf+s1k_ddi-E1_ex_0;
			s2k=s2k_chem+s2k_surf+s2k_ddi;

			//発展方程式[式(4.57)]とs場の時間発展（陽解法）の計算
			//s1ddtt=-smob1*s1k;  s1h2[i][j]=s1h[i][j]+s1ddtt*delt;
			//s2ddtt=-smob2*s2k;  s2h2[i][j]=s2h[i][j]+s2ddtt*delt;
			s1ddtt=-smob1*s1k;  s1h2[i*ND+j]=s1h[i*ND+j]+s1ddtt*delt;
			s2ddtt=-smob2*s2k;  s2h2[i*ND+j]=s2h[i*ND+j]+s2ddtt*delt;
		}
	}

	//s場を主配列に移動および変域に関する数値計算誤差の補正
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//s1h[i][j]=s1h2[i][j]; s2h[i][j]=s2h2[i][j];
			s1h[i*ND+j]=s1h2[i*ND+j]; s2h[i*ND+j]=s2h2[i*ND+j];
			//if(s1h[i][j]>=Peq){s1h[i][j]=Peq;}  if(s1h[i][j]<=-Peq){s1h[i][j]=-Peq;}
			if(s1h[i*ND+j]>=Peq){s1h[i*ND+j]=Peq;}  if(s1h[i*ND+j]<=-Peq){s1h[i*ND+j]=-Peq;}
			//if(s2h[i][j]>=Peq){s2h[i][j]=Peq;}  if(s2h[i][j]<=-Peq){s2h[i][j]=-Peq;}
			if(s2h[i*ND+j]>=Peq){s2h[i*ND+j]=Peq;}  if(s2h[i*ND+j]<=-Peq){s2h[i*ND+j]=-Peq;}
			//ss1ss2=sqrt(s1h[i][j]*s1h[i][j]+s2h[i][j]*s2h[i][j]);
			ss1ss2=sqrt(s1h[i*ND+j]*s1h[i*ND+j]+s2h[i*ND+j]*s2h[i*ND+j]);
			//if(ss1ss2>=Peq){s1h[i][j]=s1h[i][j]/ss1ss2*Peq; s2h[i][j]=s2h[i][j]/ss1ss2*Peq;}
			if(ss1ss2>=Peq){s1h[i*ND+j]=s1h[i*ND+j]/ss1ss2*Peq; s2h[i*ND+j]=s2h[i*ND+j]/ss1ss2*Peq;}
		}
	}

	//時間を進め、終了時間でない時にはスタートへ
	//if(keypress()){return 0;}//キー待ち状態
	time1=time1+1.0;								//計算カウント数の加算
	if(time1<time1max){goto start;}	//最大カウント数に到達したかどうかの判断

	end:;
  return 0;
}

//************ 初期場の設定サブル－チン *************
void ini000_s12(double *s1h, double *s2h, int ND)
{
	int i, j;	//整数
 	//srand(time(NULL)); // 乱数初期化
	int nd=ND, ndm=ND-1, nd2=ND/2;

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//s1h[i][j]=0.01*(2.0*DRND(1)-1.0);//場を最大±1%の乱数にて設定
			//s2h[i][j]=0.01*(2.0*DRND(1)-1.0);
			s1h[i*ND+j]=0.01*(2.0*DRND(1)-1.0);//場を最大±1%の乱数にて設定
			s2h[i*ND+j]=0.01*(2.0*DRND(1)-1.0);
		}
	}
}

//******* Sin, Cos のテーブルおよびビット反転テーブルの設定 ***************
void table(double *s, double *c, int *ik, int ND, int ig)
{
	int it, it1, it2, mc, mn;
	double q;
	int nd=ND, ndm=ND-1, nd2=ND/2;

	q=2.*PI/nd;
	for(it=0;it<=nd2-1;it++){
		c[it]=cos(q*it);  s[it]=sin(q*it);//Sin, Cos のテーブル
	}
	ik[0]=0;
	mn=nd2;
	mc=1;
	for(it1=1;it1<=ig;it1++){
		for(it2=0;it2<=mc-1;it2++){
			ik[it2+mc]=ik[it2]+mn;//ビット反転テーブル
		}
		mn=mn/2;
		mc=2*mc;
	}
}

//********** １次元高速フーリエ変換 **************************************
void fft(double *xrf, double *xif, double *s, double *c, int ND, int ig)
{
	int ix, ka, kb, l2, lf, mf, n2, nf;
	double tj, tr;
	int nd=ND, ndm=ND-1, nd2=ND/2;

	l2=1;
	for(lf=1;lf<=ig;lf++){
		n2=nd2/l2;
		for(mf=1;mf<=l2;mf++){
			for(nf=0;nf<=n2-1;nf++){
				ix=nf*l2;
				ka=nf+2*n2*(mf-1);
				kb=ka+n2;
				tr=xrf[ka]-xrf[kb];
				tj=xif[ka]-xif[kb];
				xrf[ka]=xrf[ka]+xrf[kb];
				xif[ka]=xif[ka]+xif[kb];
				xrf[kb]=tr*c[ix]-tj*qs*s[ix];
				xif[kb]=tj*c[ix]+tr*qs*s[ix];
			}
		}
		l2=l2*2;
	}

}

//************ ２次元高速フーリエ変換 ***********************************
void rcfft(double *xrf, double *xif, double *xr, double *xi, int *ik, double *s, double *c, int ND, int ig)
{
	int i, ic, ir, j;
	int nd=ND, ndm=ND-1, nd2=ND/2;

	for(ir=0;ir<=ndm;ir++){
		for(ic=0;ic<=ndm;ic++){
			//xrf[ic]=xr[ir][ic];  xif[ic]=xi[ir][ic];
			xrf[ic]=xr[ir*ND+ic];
			xif[ic]=xi[ir*ND+ic];
		}
		fft(xrf, xif, s, c, ND, ig);
		for(ic=0;ic<=ndm;ic++){
			//xr[ir][ic]=xrf[ik[ic]];  xi[ir][ic]=xif[ik[ic]];
			xr[ir*ND+ic]=xrf[ik[ic]];
			xi[ir*ND+ic]=xif[ik[ic]];
		}
	}
	for(ic=0;ic<=ndm;ic++){
		for(ir=0;ir<=ndm;ir++){
			//xrf[ir]=xr[ir][ic];  xif[ir]=xi[ir][ic];
			xrf[ir]=xr[ir*ND+ic];
			xif[ir]=xi[ir*ND+ic];
		}
		fft(xrf, xif, s, c, ND, ig);
		for(ir=0;ir<=ndm;ir++){
			//xr[ir][ic]=xrf[ik[ir]];  xi[ir][ic]=xif[ik[ir]];
			xr[ir*ND+ic]=xrf[ik[ir]];
			xi[ir*ND+ic]=xif[ik[ir]];
		}
	}
	if(qs>0.){return;}
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//xr[i][j]=xr[i][j]/nd/nd;  xi[i][j]=xi[i][j]/nd/nd;
			xr[i*ND+j]=xr[i*ND+j]/nd/nd;
			xi[i*ND+j]=xi[i*ND+j]/nd/nd;
		}
	}

}

//************ データ保存サブルーチン *******************************
void datsave(double *s1h, double *s2h, int ND)
{
	FILE		*stream;	//ストリームのポインタ設定
	int 		i, j;			//整数
	int nd=ND, ndm=ND-1, nd2=ND/2;

	stream = fopen("test.dat", "a");	//書き込み先のファイルを追記方式でオープン
	fprintf(stream, "%e \n", time1);	//計算カウント数の保存
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//fprintf(stream, "%e  %e  ", s1h[i][j], s2h[i][j]);	//場のデータ保存
			fprintf(stream, "%e  %e  ", s1h[i*ND+j], s2h[i*ND+j]);	//場のデータ保存
		}
	}
	fprintf(stream, "\n");	//改行の書き込み
	fclose(stream);					//ファイルをクローズ
}

void datsave_paraview(double *s1h, double *s2h, int ND)
{
	FILE	*fp;
	char	fName[256];
	int 	i, j;
	int nd=ND, ndm=ND-1, nd2=ND/2;
	
	iout = iout + 1;
	printf("dip_result%06d.vtk \n",iout);
	sprintf(fName,"dip_result%06d.vtk",iout);
	fp = fopen(fName, "w");
	fprintf(fp,"# vtk DataFile Version 3.0 \n");
	fprintf(fp,"output.vtk \n");
	fprintf(fp,"ASCII \n");
	fprintf(fp,"DATASET STRUCTURED_POINTS \n");
	fprintf(fp,"DIMENSIONS %5d %5d %5d \n",(ndm+1),(ndm+1),(1));
	fprintf(fp,"ORIGIN 0.0 0.0 0.0 \n");
	fprintf(fp,"ASPECT_RATIO 1 1 1 \n");
	fprintf(fp,"POINT_DATA %16d \n",((ndm+1)*(ndm+1)*1));
	fprintf(fp,"SCALARS Phase_field_1 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			//fprintf(stream, "%e  %e  ", s1h[i][j], s2h[i][j]);	//場のデータ保存
			//fprintf(fp,"%10.6f\n", s1h[i][j]);
			fprintf(fp,"%10.6f\n", s1h[i*ND+j]);
		}
	}
	fprintf(fp,"SCALARS Phase_field_2 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			//fprintf(stream, "%e  %e  ", s1h[i][j], s2h[i][j]);	//場のデータ保存
			//fprintf(fp,"%10.6f\n", s2h[i][j]);
			fprintf(fp,"%10.6f\n", s2h[i*ND+j]);
		}
	}
	fclose(fp);
}
//************ データ読み込みサブルーチン *******************************
void datin(double *s1h, double *s2h, int ND)
{
	FILE		*datin0;	//ストリームのポインタ設定
	int 		i, j;			//整数
	double time_ext;	//計算カウント数
	int nd=ND, ndm=ND-1, nd2=ND/2;

	datin0 = fopen("test.dat", "r");	//読み込み元のファイルをオープン
	fscanf(datin0, "%lf", &time_ext);		//計算カウント数の読込
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//fscanf(datin0, "%lf  %lf  ", &s1h[i][j], &s2h[i][j]);	//場のデータ読込
			fscanf(datin0, "%lf  %lf  ", &s1h[i*ND+j], &s2h[i*ND+j]);	//場のデータ読込
		}
	}
	fclose(datin0);//ファイルをクローズ
}