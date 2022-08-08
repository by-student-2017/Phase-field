#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

#define DRND(x) ((double)(x)/RAND_MAX*rand())//乱数の関数設定
//#define ND 256					//差分計算における計算領域一辺の分割数(高速フーリエ変換を用いるため２のべき乗)
//#define IG 8					//2^IG=ND

	//int nd=ND, ndm=ND-1; 		//計算領域の一辺の差分分割数(差分ブロック数)、ND-1を定義
	//int nd2=ND/2;				 //ND/2を定義：高速フ－リエ変換で使用
	//int ig=IG;					//2^ig=ND
	double PI=3.14159;			//円周率
	double rr=8.3145;			//ガス定数
	double time1;				//計算カウント数(時間に比例)
	int Nstep, iout;

	//double s1h[ND][ND], s2h[ND][ND];	//マルテンサイトのフェーズフィールド

	double qs;					//フ－リエ変換(qs:-1)とフ－リエ逆変換(qs:1)の区別
	//double xi[ND][ND], xr[ND][ND], xif[ND], xrf[ND];//フ－リエ変換の実部・虚部配列
	//double s[ND],c[ND];			//sinとcosのテーブル
	//int ik[ND];					//ビット反転テーブル

	void ini000(double *s1h, double *s2h, int ND);	//初期場の設定サブル－チン
	void table(double *s, double *c, int *ik, int ND, int ig);//sinとcosのテーブルとビット反転テーブルの作成サブル－チン
	void fft(double *xrf, double *xif, double *s, double *c, int ND, int ig);//１次元高速フーリエ変換
	void rcfft(double *xrf, double *xif, double *xr, double *xi, int *ik, double *s, double *c, int ND, int ig);//２次元高速フーリエ変換
	void datsave(double *s1h, double *s2h, int ND);			//デ－タ保存サブル－チン
	void datsave_paraview(double *s1h, double *s2h, int ND);//デ－タ保存サブル－チン

//******* メインプログラム ******************************************
int main(void)
{
    int ND, IG;
	int nd, ndm, nd2, ig;
	
	double s1, s2;									//マルテンサイトのフェーズフィールド
	//double ep11h0[ND][ND], ep22h0[ND][ND];			//組織内の変態歪
	//double ep11qrh0[ND][ND],	ep11qih0[ND][ND];	//拘束歪変動量のフーリエ変換
	//double ep22qrh0[ND][ND],	ep22qih0[ND][ND];	//拘束歪変動量のフーリエ変換
	double s1k_chem, s1k_str;						//ポテンシャル
	//double s1k_su[ND][ND];							//ポテンシャル
	double s2k_chem, s2k_str;						//ポテンシャル
	//double s2k_su[ND][ND];							//ポテンシャル
	double c11, c12, c44, lam0, mu0, nu0; 			//弾性定数
	double c11c, c12c, c44c;			 			//弾性定数
	double eta_s1[4][4], eta_s2[4][4];				//アイゲン歪成分
	//double ec11[ND][ND], ec22[ND][ND];				//拘束歪変動量（実空間）
	double ep11T, ep22T;
	double ep11_0, ep22_0;							//組織内の変態歪の平均値
	double ep11_a, ep22_a, ep12_a, ep21_a;			//外力に起因する歪
	double sig11_a, sig22_a;						//外力
	double Z11ep, Z12ep, Z21ep, Z22ep;				//フーリエ逆変換時の係数
	double sum11, sum22;							//s1とs2の空間積分
	double s1ddtt, s2ddtt;							//s1とs2の時間変化量（発展方程式の左辺）
	double el_fac;									//弾性定数の規格化定数

	int   i, j, k, l, ii, jj, kk, iii, jjj;			//整数
	int   p, q, m, n;								//整数
	int   ip, im, jp, jm, Nstep;					//整数
	double al, temp, delt;							//計算領域、温度、時間きざみ
	double time1max;								//計算カウント数の最大値（計算終了カウント）
	double b1, vm0, atom_n;							//差分プロック１辺の長さ、モル体積、単位胞内の原子数
	double smob;									//モビリティー（結晶変態の緩和係数）
	double nxx, nyy, nxy, alnn;						//フーリエ空間の基本ベクトルの積、ノルム

	double AA0, AA1, AA2, AA3;						//ギブスエネルギー内の係数
	double AA0e;									//ギブスエネルギー内の係数
	double a1_c, b1_c, c1_c;						//格子定数
	double a1_t, b1_t, c1_t;						//格子定数
	double kappa_s1, kappa_s2;						//勾配エネルギ－係数
	double kappa_s1c, kappa_s2c;					//勾配エネルギ－係数
	double ds_fac;									//結晶変態の揺らぎ係数
	double eta1, eta2;

//****** 計算条件および物質定数の設定 ****************************************
	printf("---------------------------------\n");
	printf("read parameters form parameters.txt\n");
	FILE *fp;
	char name[40], comment[72];
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
	temp    = data[2];
	al      = data[3];	// [nm]
	smob    = data[4];
	ds_fac  = data[5];
	AA0e    = data[6];
	AA1     = data[7];
	AA2     = data[8];
	AA3     = data[9];
	kappa_s1c=data[10];
	kappa_s2c=data[11];
	vm0     = data[12];
	time1max= int(data[13]);
	Nstep   = int(data[14]);
	c11c    = data[15];
	c12c    = data[16];
	c44c    = data[17];
	lam0    = data[18];
	mu0     = data[19];
	eta1    = data[20];
	eta2    = data[21];
	sig22_a = data[22];
	printf("---------------------------------\n");
	//
	IG=int(log2(ND));
	nd=ND; 					
	ndm=ND-1;				//define ND-1
	nd2=ND/2;				//define ND/2 for FFT
	ig=IG;					//2^ig=ND
	//
	double *s1h      = (double *)malloc(sizeof(double)*( ND*ND + ND ));//マルテンサイトのフェーズフィールド
	double *s2h      = (double *)malloc(sizeof(double)*( ND*ND + ND ));//マルテンサイトのフェーズフィールド
	//
	double *ep11h0   = (double *)malloc(sizeof(double)*( ND*ND + ND ));//組織内の変態歪
	double *ep22h0   = (double *)malloc(sizeof(double)*( ND*ND + ND ));//組織内の変態歪
	double *ep11qrh0 = (double *)malloc(sizeof(double)*( ND*ND + ND ));//拘束歪変動量のフーリエ変換
	double *ep11qih0 = (double *)malloc(sizeof(double)*( ND*ND + ND ));//拘束歪変動量のフーリエ変換
	double *ep22qrh0 = (double *)malloc(sizeof(double)*( ND*ND + ND ));//拘束歪変動量のフーリエ変換
	double *ep22qih0 = (double *)malloc(sizeof(double)*( ND*ND + ND ));//拘束歪変動量のフーリエ変換
	//
	double *s1k_su   = (double *)malloc(sizeof(double)*( ND*ND + ND ));//ポテンシャル
	double *s2k_su   = (double *)malloc(sizeof(double)*( ND*ND + ND ));//ポテンシャル
	//
	double *ec11     = (double *)malloc(sizeof(double)*( ND*ND + ND ));//拘束歪変動量（実空間）
	double *ec22     = (double *)malloc(sizeof(double)*( ND*ND + ND ));//拘束歪変動量（実空間）
	//
	double *xi       = (double *)malloc(sizeof(double)*( ND*ND + ND ));//フ－リエ変換の実部・虚部配列
	double *xr       = (double *)malloc(sizeof(double)*( ND*ND + ND ));//フ－リエ変換の実部・虚部配列
	//
	double *xif      = (double *)malloc(sizeof(double)*( ND ));//フ－リエ変換の実部・虚部配列
	double *xrf      = (double *)malloc(sizeof(double)*( ND ));//フ－リエ変換の実部・虚部配列
	//
	double *s        = (double *)malloc(sizeof(double)*( ND ));//sinのテーブル
	double *c        = (double *)malloc(sizeof(double)*( ND ));//cosのテーブル
	//
	int *ik       = (int *)malloc(sizeof(int)*( ND ));//ビット反転テーブル
	//
	//printf("DELT(0.2)=  ");	scanf(" %lf",&delt);	//時間きざみ入力	//delt=0.2;

	//temp=500.0;					//温度(K)
	//al=500.0*1.0E-09;			//計算領域(m)
	al=al*1.0E-09;				//計算領域(m)
	b1=al/nd;					//差分ブロックの長さ

	time1=0.0;					//初期計算カウント数の設定
	//time1max=1.0+1.0e+07;		//最大計算カウント数の設定

	//smob=1.0;					//モビリティー（結晶変態の緩和係数）
	//ds_fac=0.01;				//結晶変態の揺らぎ係数

	//AA0=1000.0/rr/temp;			//マルテンサイト変態の化学的駆動力
	AA0=AA0e/rr/temp;			//マルテンサイト変態の化学的駆動力
	//AA1=1.0;  AA2=3.0*AA1+12.0;  AA3=2.0*AA1+12.0;	//ギブスエネルギー内の係数

	//kappa_s1=kappa_s2=5.0e-15/rr/temp/b1/b1;	//勾配エネルギ－係数
	kappa_s1=kappa_s1c/rr/temp/b1/b1;	//勾配エネルギ－係数
	kappa_s2=kappa_s2c/rr/temp/b1/b1;	//勾配エネルギ－係数

	//a1_c=b1_c=c1_c=3.5E-10;						//格子定数(m)
	//atom_n=4.0;  vm0=6.02E23*a1_c*b1_c*c1_c/atom_n;	//モル体積の計算（fccを仮定）

//*** s1場のアイゲン歪の設定 ***************
	//eta_s1[1][1]=0.05; eta_s1[2][2]=-0.05;
	eta_s1[1][1]=eta1; eta_s1[2][2]=eta2;
	eta_s1[3][3]=0.;
	eta_s1[1][2]=eta_s1[2][1]=eta_s1[1][3]=eta_s1[3][1]=eta_s1[2][3]=eta_s1[3][2]=0.;

//*** s2場のアイゲン歪の設定 ***************
	eta_s2[1][1]=eta_s1[2][2];
	eta_s2[2][2]=eta_s1[1][1];
	eta_s2[3][3]=0.;
	eta_s2[1][2]=eta_s2[2][1]=eta_s2[1][3]=eta_s2[3][1]=eta_s2[2][3]=eta_s2[3][2]=0.;

//***** Niの弾性定数 ****************************
	el_fac=1.0E+11*vm0/rr/temp;
	//c11=2.508*el_fac;
	//c44=1.235*el_fac;
	//c12=1.500*el_fac;
	c11=c11c*el_fac;
	c12=c12c*el_fac;
	c44=c44c*el_fac;
	//c12=c11-2.0*c44;
	if(lam0==0.0){
		lam0=c12;//cijのデータからラーメ定数を設定
	}
	if(mu0==0.0){
		mu0=c44;//cijのデータからラーメ定数を設定
	}
	nu0=lam0/2.0/(lam0+mu0);//ポアソン比
	printf("nu0= %f  \n", nu0);//ポアソン比の値を表示（ラーメ定数の妥当性の確認のため）

//*** 外力の設定 *******************************
 	//sig22_a=0.0;//本計算では、外力を考慮していないので０を設定
	ep11_a=-0.5*lam0/mu0/(3.0*lam0+2.0*mu0)*sig22_a;
	ep22_a=(lam0+mu0)/mu0/(3.0*lam0+2.0*mu0)*sig22_a;
	ep12_a=ep21_a=0.0;

//*** sinおよびcosテ－ブル、ビット反転テーブル、および初期場の設定 ***************

	table(s, c, ik, ND, ig);	//sinおよびcosテ－ブルとビット反転テーブルの設定

	ini000(s1h, s2h, ND);		//初期場の設定

//**** シミュレーションスタート ******************************
//Nstep = 10;
iout = -1;
start: ;

	//if(time1<=100.){Nstep=10;} else{Nstep=200;}		//データ保存する時間間隔の変更
	//if((((int)(time1) % Nstep)==0)) {datsave(s1h, s2h, ND);} 	//一定繰返しカウント毎に組織データを保存
	if((((int)(time1) % Nstep)==0)) {datsave_paraview(s1h, s2h, ND);} 	//一定繰返しカウント毎に組織データを保存

//***** 勾配ポテンシャル ***********************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ip=i+1; im=i-1; jp=j+1; jm=j-1;
			if(i==ndm){ip=0;}  if(i==0){im=ndm;}	//周期的境界条件
			if(j==ndm){jp=0;}  if(j==0){jm=ndm;}
			//s1k_su[i][j]=-kappa_s1*(s1h[ip][j]+s1h[im][j]+s1h[i][jp]+s1h[i][jm]-4.0*s1h[i][j]);//式(4.2.4)
			//s2k_su[i][j]=-kappa_s2*(s2h[ip][j]+s2h[im][j]+s2h[i][jp]+s2h[i][jm]-4.0*s2h[i][j]);
			s1k_su[i*ND+j]=-kappa_s1*(s1h[ip*ND+j]+s1h[im*ND+j]+s1h[i*ND+jp]+s1h[i*ND+jm]-4.0*s1h[i*ND+j]);//式(4.2.4)
			s2k_su[i*ND+j]=-kappa_s2*(s2h[ip*ND+j]+s2h[im*ND+j]+s2h[i*ND+jp]+s2h[i*ND+jm]-4.0*s2h[i*ND+j]);
		}
	}

//**** アイゲン歪場[式(4.7)]のフ－リエ変換 ep11 ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//xr[i][j]=ep11h0[i][j]=eta_s1[1][1]*s1h[i][j]+eta_s2[1][1]*s2h[i][j];
			//xi[i][j]=0.0;
			xr[i*ND+j]=ep11h0[i*ND+j]=eta_s1[1][1]*s1h[i*ND+j]+eta_s2[1][1]*s2h[i*ND+j];
			xi[i*ND+j]=0.0;
		}
	}
	qs=-1.; rcfft(xrf, xif, xr, xi, ik, s, c, ND, ig);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){ 
			//ep11qrh0[i][j]=xr[i][j];
			//ep11qih0[i][j]=xi[i][j];
			ep11qrh0[i*ND+j]=xr[i*ND+j];
			ep11qih0[i*ND+j]=xi[i*ND+j];
		}
	}
	//ep11qrh0[0][0]=ep11qih0[0][0]=0.;
	ep11qrh0[0]=ep11qih0[0]=0.0;

//**** アイゲン歪場[式(4.7)]のフ－リエ変換 ep22 ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//xr[i][j]=ep22h0[i][j]=eta_s1[2][2]*s1h[i][j]+eta_s2[2][2]*s2h[i][j];
			//xi[i][j]=0.0;
			xr[i*ND+j]=ep22h0[i*ND+j]=eta_s1[2][2]*s1h[i*ND+j]+eta_s2[2][2]*s2h[i*ND+j];
			xi[i*ND+j]=0.0;
		}
	}
	qs=-1.; rcfft(xrf, xif, xr, xi, ik, s, c, ND, ig);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){ 
			//ep22qrh0[i][j]=xr[i][j];
			//ep22qih0[i][j]=xi[i][j];
			ep22qrh0[i*ND+j]=xr[i*ND+j];
			ep22qih0[i*ND+j]=xi[i*ND+j];
		}
	}
	//ep22qrh0[0][0]=ep22qih0[0][0]=0.;
	ep22qrh0[0]=ep22qih0[0]=0.0;

//*** アイゲン歪場の平均値の算出 ***
	sum11=sum22=0.;
	for(i=0;i<=ndm;i++){
		//for(j=0;j<=ndm;j++){ sum11+=ep11h0[i][j];  sum22+=ep22h0[i][j]; }
		for(j=0;j<=ndm;j++){ sum11+=ep11h0[i*ND+j];  sum22+=ep22h0[i*ND+j]; }
	}
	ep11_0=sum11/nd/nd;  ep22_0=sum22/nd/nd;

//***** 拘束歪変動量ec11の計算[式(4.9)] *************************************
	for(i=0;i<=ndm;i++){
		if(i<=nd2-1){ii=i;}  if(i>=nd2){ii=i-nd;}
		for(j=0;j<=ndm;j++){
			if(j<=nd2-1){jj=j;}  if(j>=nd2){jj=j-nd;}
			alnn=sqrt((double)ii*(double)ii+(double)jj*(double)jj);
			if(alnn==0.){alnn=1.;}
			nxx=(double)ii/alnn*(double)ii/alnn;
			nyy=(double)jj/alnn*(double)jj/alnn;
			Z11ep=nxx*(2.0*(1.0-nu0)-nxx-nu0/(1.0-nu0)*nyy)/(1.0-2.0*nu0);
			Z12ep=nxx*(2.0*nu0      -nyy-nu0/(1.0-nu0)*nxx)/(1.0-2.0*nu0);
			//xr[i][j]=Z11ep*ep11qrh0[i][j]+Z12ep*ep22qrh0[i][j];
			//xi[i][j]=Z11ep*ep11qih0[i][j]+Z12ep*ep22qih0[i][j];
			xr[i*ND+j]=Z11ep*ep11qrh0[i*ND+j]+Z12ep*ep22qrh0[i*ND+j];
			xi[i*ND+j]=Z11ep*ep11qih0[i*ND+j]+Z12ep*ep22qih0[i*ND+j];
		}
	}
	qs=1.; rcfft(xrf, xif, xr, xi, ik, s, c, ND, ig);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//ec11[i][j]=xr[i][j];
			ec11[i*ND+j]=xr[i*ND+j];
		}
	}

//***** 拘束歪変動量ec22の計算[式(4.9)] *****************************
	for(i=0;i<=ndm;i++){
		if(i<=nd2-1){ii=i;}  if(i>=nd2){ii=i-nd;}
		for(j=0;j<=ndm;j++){
			if(j<=nd2-1){jj=j;}  if(j>=nd2){jj=j-nd;}
			alnn=sqrt((double)ii*(double)ii+(double)jj*(double)jj);
			if(alnn==0.){alnn=1.;}
			nxx=(double)ii/alnn*(double)ii/alnn;
			nyy=(double)jj/alnn*(double)jj/alnn;
			Z21ep=nyy*(2.0*nu0      -nxx-nu0/(1.0-nu0)*nyy)/(1.0-2.0*nu0);
			Z22ep=nyy*(2.0*(1.0-nu0)-nyy-nu0/(1.0-nu0)*nxx)/(1.0-2.0*nu0);
			//xr[i][j]=Z21ep*ep11qrh0[i][j]+Z22ep*ep22qrh0[i][j];
			//xi[i][j]=Z21ep*ep11qih0[i][j]+Z22ep*ep22qih0[i][j];
			xr[i*ND+j]=Z21ep*ep11qrh0[i*ND+j]+Z22ep*ep22qrh0[i*ND+j];
			xi[i*ND+j]=Z21ep*ep11qih0[i*ND+j]+Z22ep*ep22qih0[i*ND+j];
		}
	}
	qs=1.; rcfft(xrf, xif, xr, xi, ik, s, c, ND, ig);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//ec22[i][j]=xr[i][j];
			ec22[i*ND+j]=xr[i*ND+j];
		}
	}

//******  ポテンシャルの計算 ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){

			//s1=s1h[i][j];  	s2=s2h[i][j];
			s1=s1h[i*ND+j];  	s2=s2h[i*ND+j];

//******  化学ポテンシャルの計算[式(4.4)] ********************************
			s1k_chem=AA0*s1*(AA1-AA2*s1+AA3*(s1*s1+s2*s2));
			s2k_chem=AA0*s2*(AA1-AA2*s2+AA3*(s1*s1+s2*s2));

//******  弾性ポテンシャルの計算[式(4.8)] ********************************

			//ep11T=ep11h0[i][j]-ep11_0-ec11[i][j]-ep11_a;
			//ep22T=ep22h0[i][j]-ep22_0-ec22[i][j]-ep22_a;
			ep11T=ep11h0[i*ND+j]-ep11_0-ec11[i*ND+j]-ep11_a;
			ep22T=ep22h0[i*ND+j]-ep22_0-ec22[i*ND+j]-ep22_a;

			s1k_str=ep11T*((lam0+2.0*mu0)*eta_s1[1][1]+lam0*eta_s1[2][2])
				   +ep22T*((lam0+2.0*mu0)*eta_s1[2][2]+lam0*eta_s1[1][1]);
			s2k_str=ep11T*((lam0+2.0*mu0)*eta_s2[1][1]+lam0*eta_s2[2][2])
				   +ep22T*((lam0+2.0*mu0)*eta_s2[2][2]+lam0*eta_s2[1][1]);

//****** フェーズフィールドの時間発展の計算[式(4.10)] ********************************
			//s1ddtt=-smob*(s1k_chem+s1k_su[i][j]+s1k_str);
			//s2ddtt=-smob*(s2k_chem+s2k_su[i][j]+s2k_str);
			//s1h[i][j]=s1h[i][j]+( s1ddtt+ds_fac*(2.0*DRND(1.)-1.0) )*delt;//陽解法
			//s2h[i][j]=s2h[i][j]+( s2ddtt+ds_fac*(2.0*DRND(1.)-1.0) )*delt;
			s1ddtt=-smob*(s1k_chem+s1k_su[i*ND+j]+s1k_str);
			s2ddtt=-smob*(s2k_chem+s2k_su[i*ND+j]+s2k_str);
			s1h[i*ND+j]=s1h[i*ND+j]+( s1ddtt+ds_fac*(2.0*DRND(1.)-1.0) )*delt;//陽解法
			s2h[i*ND+j]=s2h[i*ND+j]+( s2ddtt+ds_fac*(2.0*DRND(1.)-1.0) )*delt;

//*** sの変域(0<=s<=1)の補正 ***
			//if(s1h[i][j]>=1.0){s1h[i][j]=1.0;}  if(s1h[i][j]<=0.0){s1h[i][j]=0.0;}
			//if(s2h[i][j]>=1.0){s2h[i][j]=1.0;}  if(s2h[i][j]<=0.0){s2h[i][j]=0.0;}
			if(s1h[i*ND+j]>=1.0){s1h[i*ND+j]=1.0;}
			if(s1h[i*ND+j]<=0.0){s1h[i*ND+j]=0.0;}
			if(s2h[i*ND+j]>=1.0){s2h[i*ND+j]=1.0;}
			if(s2h[i*ND+j]<=0.0){s2h[i*ND+j]=0.0;}
		}
	}

	//if(keypress()){return 0;}//キー待ち状態

	time1=time1+1.0;								//計算カウント数の加算
	if(time1<time1max){goto start;}	//最大カウント数に到達したかどうかの判断

end:;
  return 0;
}

//************ 初期場の設定サブル－チン *************
void ini000(double *s1h, double *s2h, int ND)
{
	int i, j;
	//srand(time(NULL)); // 乱数初期化
	int nd=ND, ndm=ND-1, nd2=ND/2;

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//s1h[i][j]=0.2*DRND(1.0); s2h[i][j]=0.2*DRND(1.0);//場を最大20%の乱数にて設定
			//if(abs(j-nd2)<(nd/40)){s1h[i][j]=DRND(1.0); s2h[i][j]=DRND(1.0);}
			s1h[i*ND+j]=0.2*DRND(1.0);//場を最大20%の乱数にて設定
			s2h[i*ND+j]=0.2*DRND(1.0);//場を最大20%の乱数にて設定
		}
	}
}

//******* Sin, Cos のテーブルおよびビット反転テーブルの設定 ***************
void table(double *s, double *c, int *ik, int ND, int ig)
{
	int it, it1, it2, mc, mn;
	double q;
	int nd=ND, ndm=ND-1, nd2=ND/2;

	q=2.0*PI/nd;
	for(it=0;it<=nd2-1;it++){ c[it]=cos(q*it); s[it]=sin(q*it); }//Sin, Cos のテーブル

	ik[0]=0; mn=nd2; mc=1;
	for(it1=1;it1<=ig;it1++){
		for(it2=0;it2<=mc-1;it2++){
			ik[it2+mc]=ik[it2]+mn;				//ビット反転テーブル
		}
		mn=mn/2; mc=2*mc;
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
				tr=xrf[ka]-xrf[kb];  					tj=xif[ka]-xif[kb];
				xrf[ka]=xrf[ka]+xrf[kb]; 			xif[ka]=xif[ka]+xif[kb];
				xrf[kb]=tr*c[ix]-tj*qs*s[ix];	xif[kb]=tj*c[ix]+tr*qs*s[ix];
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
			//xrf[ic]=xr[ir][ic];	xif[ic]=xi[ir][ic];
			xrf[ic]=xr[ir*ND+ic];
			xif[ic]=xi[ir*ND+ic];
		}
		fft(xrf, xif, s, c, ND, ig);
		for(ic=0;ic<=ndm;ic++){
			//xr[ir][ic]=xrf[ik[ic]];	xi[ir][ic]=xif[ik[ic]];
			xr[ir*ND+ic]=xrf[ik[ic]];
			xi[ir*ND+ic]=xif[ik[ic]];
		}
	}
	for(ic=0;ic<=ndm;ic++){
		for(ir=0;ir<=ndm;ir++){
			//xrf[ir]=xr[ir][ic];	xif[ir]=xi[ir][ic];
			xrf[ir]=xr[ir*ND+ic];
			xif[ir]=xi[ir*ND+ic];
		}
		fft(xrf, xif, s, c, ND, ig);
		for(ir=0;ir<=ndm;ir++){
			//xr[ir][ic]=xrf[ik[ir]];	xi[ir][ic]=xif[ik[ir]];
			xr[ir*ND+ic]=xrf[ik[ir]];
			xi[ir*ND+ic]=xif[ik[ir]];
		}
	}
	if(qs>0.0){return;}
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//xr[i][j]=xr[i][j]/nd/nd;	xi[i][j]=xi[i][j]/nd/nd;
			xr[i*ND+j]=xr[i*ND+j]/nd/nd;
			xi[i*ND+j]=xi[i*ND+j]/nd/nd;
		}
	}

}

//************ データ保存サブルーチン *******************************
void datsave(double *s1h, double *s2h, int ND)
{
	FILE *stream;		//ストリームのポインタ設定
	int i, j;			//整数
	int nd=ND, ndm=ND-1, nd2=ND/2;

	stream = fopen("test.dat", "a");	//書き込み先のファイルを追記方式でオープン
	fprintf(stream, "%e\n", time1);		//計算カウント数の保存
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//fprintf(stream, "%e  %e ", s1h[i][j], s2h[i][j]);//場のデータ保存
			fprintf(stream, "%e  %e ", s1h[i*ND+j], s2h[i*ND+j]);
		}
	}
	fprintf(stream, "\n");	//改行の書き込み
	fclose(stream);			//ファイルをクローズ
}

void datsave_paraview(double *s1h, double *s2h, int ND)
{
	FILE	*fp;
	char	fName[256];
	int 	i, j;
	int nd=ND, ndm=ND-1, nd2=ND/2;
	
	iout = iout + 1;
	printf("pf_result%06d.vtk \n",iout);
	sprintf(fName,"pf_result%06d.vtk",iout);
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
			//fprintf(stream, "%e  %e ", s1h[i][j], s2h[i][j]);//場のデータ保存
			//fprintf(fp,"%10.6f\n", s1h[i][j]);
			fprintf(fp,"%10.6f\n", s1h[i*ND+j]);
		}
	}
	fprintf(fp,"SCALARS Phase_field_2 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			//fprintf(stream, "%e  %e ", s1h[i][j], s2h[i][j]);//場のデータ保存
			//fprintf(fp,"%10.6f\n", s2h[i][j]);
			fprintf(fp,"%10.6f\n", s2h[i*ND+j]);
		}
	}
	fclose(fp);
}