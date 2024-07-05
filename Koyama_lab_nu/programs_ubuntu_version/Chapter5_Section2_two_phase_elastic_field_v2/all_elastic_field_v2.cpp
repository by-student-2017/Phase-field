//フェーズフィールドのデータ(ph.dat)を読んで、
//弾性場を計算するプログラム。
//各弾性場は、一括してel_field.datに保存される。
//同時に、各弾性場の２次元画像も保存される。

#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

#define DRND(x) ((double)(x)/RAND_MAX*rand())//乱数の関数設定
//#define ND 128					//差分計算における計算領域一辺の分割数(高速フーリエ変換を用いるため２のべき乗)
//#define IG 7					//2^IG=ND

	//int nd=ND, ndm=ND-1; 		//計算領域の一辺の差分分割数(差分ブロック数)、ND-1を定義
	//int nd2=ND/2;				//ND/2を定義：高速フ－リエ変換で使用
	//int ig=IG;					//2^ig=ND
	double PI=3.14159;			//円周率
	double rr=8.3145;			//ガス定数
	//double time1;				//計算カウント数(時間に比例)
	int iout=-1;

	double c0;					//フェーズフィールドの平均値
	//double ch[ND][ND];			//場（フェーズフィールド）
	double eta_c[4][4];			//格子ミスマッチ
	double cec[4][4][4][4];		//弾性定数
	double sigma[4][4];			//アイゲン応力（アイゲンひずみ分だけ弾性変形した時の応力）
	//double ep11c[ND][ND], ep22c[ND][ND], ep12c[ND][ND];									//全ひずみ
	//double sig11[ND][ND], sig22[ND][ND], sig33[ND][ND], sig12[ND][ND];	//弾性応力
	//double Estr[ND][ND];		//弾性ひずみエネルギー
	//double u1[ND][ND], u2[ND][ND], u3[ND][ND];	//変位
	//double fld1[ND][ND];		//場の受け渡し用の作業行列

	int qs;					//フ－リエ変換(qs:-1)とフ－リエ逆変換(qs:1)の区別
	//double xi[ND][ND], xr[ND][ND], xif[ND], xrf[ND];//フ－リエ変換の実部・虚部配列
	//double s[ND],c[ND];	//sinとcosのテーブル
	//int ik[ND];					//ビット反転テーブル

	void table(double *s, double *c, int *ik, int ND, int ig);	//sinとcosのテーブルとビット反転テーブルの作成サブル－チン
	void fft(double *xrf, double *xif, double *s, double *c, int ND, int ig);	//１次元高速フーリエ変換
	void rcfft(double *xrf, double *xif, double *xr, double *xi, int *ik, double *s, double *c, int ND, int ig);	//２次元高速フーリエ変換
	double zcij(int i0, int j0, int iii, int jjj, int ND);//弾性エネルギー関数の係数計算（フーリエ空間）
	double zuij(int i0, int j0, int iii, int ND);//変位場の係数計算（フーリエ空間）

	void datin(double *ch, int ND);	//初期場読み込み用サブル－チン
	void datsave(double *ch, double *Estr, 
			 double *ep11c, double *ep22c, double *ep12c,
			 double *sig11, double *sig22, double *sig33, double *sig12, 
			 double *u1, double *u2, int ND);	//デ－タ保存サブル－チン
	void datsave_paraview(double *ch, double *Estr, 
			 double *ep11c, double *ep22c, double *ep12c,
			 double *sig11, double *sig22, double *sig33, double *sig12, 
			 double *u1, double *u2, int ND);	//デ－タ保存サブル－チン

//******* メインプログラム ****************************************************
int main(void)
{
    int ND, IG;
	int ndm, ig;
	
	//int i, j, k, l, ii, jj, kk, iii, jjj, ief;			//整数
	int i, j, k, l, iii, jjj;	//integer
	//int ip, im, jp, jm;									//整数

	//double c;	//場（フェーズフィールド）
	//double qrh1[ND][ND], qih1[ND][ND];					//場のフーリエ変換
	//double ec11[ND][ND], ec22[ND][ND], ec12[ND][ND];	//拘束歪配列

	//double al, temp, delt;								//計算領域の１辺の長さ、温度、時間きざみ
	//double time1max;									//計算カウント数の最大値（計算を止める際に使用）
	//double b1, vm0, atom_n;								//差分ブロックの１辺の長さ、モル体積、単位胞内の原子数
	//double a1;											//Feの格子定数

	double el_mag;										//弾性率に関する作業変数
	double c11, c44, c12; 								//弾性率
	double ep000;										//格子ミスマッチ
	//double epc11, epc22, epc33, epc12, epc13, epc23;	//拘束ひずみ
	//double ep011, ep022, ep033, ep012, ep013, ep023;	//アイゲンひずみ
	//double dep011, dep022, dep033, dep012, dep023, dep013;
	//double epc011, epc022, epc033, epc012, epc023, epc013;
	//double ef1, ef2, ef11, ef22, ef12;					//弾性関数

	//double Estr1, Estr2, Estr3, Estr4, Estr5;			//弾性ひずみエネルギー計算の際の作業変数
	double Estr1, Estr2, Estr3;							//弾性ひずみエネルギー計算の際の作業変数

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
	// pfm1 series
	ND      = int(data[0]);
	//c2a     = data[1];
	//delt    = data[2];	// timestep
	//temp    = data[3];	// [K]
	//al      = data[4];	// [nm]
	//cmob22  = data[5];
	//om_12e  = data[6];
	//kapa_c2c= data[7];
	//time1max= data[8];
	//Nstep   = int(data[9]);
	ep000   = data[10];
	el_mag  = data[11];
	c11     = data[12];
	c12     = data[13];
	c44     = data[14];
	printf("---------------------------------\n");
	//
	IG=int(log2(ND));
	ndm=ND-1;				//define ND-1
	ig=IG;					//2^ig=ND
	//
	//double ch[(nstep+1)*nstep]; // [(nstep+1)*nstep]=[nstep*nstep+nstep], a[i][j]= a[i*n+j] for a[][n]
	//a[z * ySize * xSize + y * xSize + x], a[i][j][k]=a[i*n*m + j*n + k]
	double *ch    = (double *)malloc(sizeof(double)*( ND*ND ));	//phase field
	//
	double *ep11c = (double *)malloc(sizeof(double)*( ND*ND ));	//Strain
	double *ep22c = (double *)malloc(sizeof(double)*( ND*ND ));
	//double *ep33c = (double *)malloc(sizeof(double)*( ND*ND ));
	double *ep12c = (double *)malloc(sizeof(double)*( ND*ND ));
	//double *ep13c = (double *)malloc(sizeof(double)*( ND*ND ));
	//double *ep23c = (double *)malloc(sizeof(double)*( ND*ND ));
	//
	double *sig11 = (double *)malloc(sizeof(double)*( ND*ND ));	//Elastic stress
	double *sig22 = (double *)malloc(sizeof(double)*( ND*ND ));
	double *sig33 = (double *)malloc(sizeof(double)*( ND*ND ));
	double *sig12 = (double *)malloc(sizeof(double)*( ND*ND ));
	//double *sig13 = (double *)malloc(sizeof(double)*( ND*ND ));
	//double *sig23 = (double *)malloc(sizeof(double)*( ND*ND ));
	//
	double *Estr  = (double *)malloc(sizeof(double)*( ND*ND ));	//Elastic strain energy density
	//
	double *u1    = (double *)malloc(sizeof(double)*( ND*ND ));	//Displacement
	double *u2    = (double *)malloc(sizeof(double)*( ND*ND ));
	//double *u3    = (double *)malloc(sizeof(double)*( ND*ND ));
	//
	double *qrh1  = (double *)malloc(sizeof(double)*( ND*ND ));	//Fourier transform of the field
	double *qih1  = (double *)malloc(sizeof(double)*( ND*ND ));
	//
	double *ec11  = (double *)malloc(sizeof(double)*( ND*ND ));	//Constrained strain array
	double *ec22  = (double *)malloc(sizeof(double)*( ND*ND ));
	//double *ec33  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *ec12  = (double *)malloc(sizeof(double)*( ND*ND ));
	//double *ec13  = (double *)malloc(sizeof(double)*( ND*ND ));
	//double *ec23  = (double *)malloc(sizeof(double)*( ND*ND ));
	//
	double *xr    = (double *)malloc(sizeof(double)*( ND*ND ));	//Real or imaginary array of Fourier transform
	double *xi    = (double *)malloc(sizeof(double)*( ND*ND ));
	double *xrf   = (double *)malloc(sizeof(double)*( ND ));
	double *xif   = (double *)malloc(sizeof(double)*( ND ));
	//
	double *s     = (double *)malloc(sizeof(double)*( ND ));	//sin and cos table
	double *c     = (double *)malloc(sizeof(double)*( ND ));
	//
	int *ik       = (int *)malloc(sizeof(int)*( ND ));	//Bit inversion table for FFT
	
	//
	//temp=900.;				//温度(K)
	//al=100.*1.0E-09;	//計算領域(m)
	//b1=al/nd;					//差分ブロックの長さ

	//time1=0.0;						//初期計算カウント数の設定
	//time1max=1.0+1.0e+07;	//最大計算カウント数の設定

	//a1=2.8664E-10; 	//Fe(bcc)の格子定数（実験デ－タ）
	//atom_n=2.;	vm0=6.02E23*a1*a1*a1/atom_n;	//モル体積の計算（bccをの場合）

//*** 格子ミスマッチの設定 ***
	//ep000=0.05;
	eta_c[1][1]=eta_c[2][2]=eta_c[3][3]=ep000; 			//格子ミスマッチ
	eta_c[1][2]=eta_c[2][1]=eta_c[1][3]=eta_c[3][1]=eta_c[2][3]=eta_c[3][2]=0.0;

//***** Fe(bcc)の弾性定数 ****************************
  //el_mag=1.0E+11/1.0E+09;
  c11=2.33*el_mag;
  c12=1.35*el_mag;
  c44=1.18*el_mag;

	for(i=0;i<=3;i++){
		for(j=0;j<=3;j++){
			for(k=0;k<=3;k++){
				for(l=0;l<=3;l++){
					cec[i][j][k][l]=0.0;//弾性定数配列の初期化
				}
			}
		}
	}

	cec[1][1][1][1]=c11;
	cec[2][2][2][2]=c11;
	cec[3][3][3][3]=c11;
	cec[1][2][1][2]=cec[1][2][2][1]=cec[2][1][1][2]=cec[2][1][2][1]=c44;
	cec[2][3][2][3]=cec[2][3][3][2]=cec[3][2][2][3]=cec[3][2][3][2]=c44;
	cec[1][3][1][3]=cec[1][3][3][1]=cec[3][1][1][3]=cec[3][1][3][1]=c44;
	cec[1][1][2][2]=cec[2][2][1][1]=c12;
	cec[1][1][3][3]=cec[3][3][1][1]=c12;
	cec[2][2][3][3]=cec[3][3][2][2]=c12;

//--- アイゲン応力（アイゲンひずみ分だけ弾性変形した時の応力）--------------
	sigma[1][1]=cec[1][1][1][1]*eta_c[1][1]
			   +cec[1][1][2][2]*eta_c[2][2]
			   +cec[1][1][3][3]*eta_c[3][3]; //sigma1
	sigma[2][2]=cec[2][2][1][1]*eta_c[1][1]
			   +cec[2][2][2][2]*eta_c[2][2]
			   +cec[2][2][3][3]*eta_c[3][3]; //sigma2
	sigma[3][3]=cec[3][3][1][1]*eta_c[1][1]
			   +cec[3][3][2][2]*eta_c[2][2]
			   +cec[3][3][3][3]*eta_c[3][3]; //sigma3
	sigma[2][3]=cec[2][3][2][3]*eta_c[2][3]+cec[2][3][3][2]*eta_c[3][2]; //sigma4
	sigma[3][2]=cec[3][2][2][3]*eta_c[2][3]+cec[3][2][3][2]*eta_c[3][2]; //sigma4
	sigma[1][3]=cec[1][3][1][3]*eta_c[1][3]+cec[1][3][3][1]*eta_c[3][1]; //sigma5
	sigma[3][1]=cec[3][1][1][3]*eta_c[1][3]+cec[3][1][3][1]*eta_c[3][1]; //sigma5
	sigma[1][2]=cec[1][2][1][2]*eta_c[1][2]+cec[1][2][2][1]*eta_c[2][1]; //sigma6
	sigma[2][1]=cec[2][1][1][2]*eta_c[1][2]+cec[2][1][2][1]*eta_c[2][1]; //sigma6

//*** sinおよびcosテ－ブル、ビット反転テーブル、および初期場の設定 ***************

	table(s, c, ik, ND, ig);
	//ini000();
	datin(ch, ND); //フェーズフィールドの読込

//**** 弾性場解析の計算スタート ******************************
//start: ;

//**** 場（フェーズフィールド）のフ－リエ変換 ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//xr[i][j]=ch[i][j]; xi[i][j]=0.;
			xr[i*ND+j]=ch[i*ND+j];
			xi[i*ND+j]=0.0;
		}
	}
	qs=-1; rcfft(xrf, xif, xr, xi, ik, s, c, ND, ig);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//qrh1[i][j]=xr[i][j]; qih1[i][j]=xi[i][j];
			qrh1[i*ND+j]=xr[i*ND+j];
			qih1[i*ND+j]=xi[i*ND+j];
		}
	}
	//qrh1[0][0]=qih1[0][0]=0.;
	qrh1[0]=qih1[0]=0.0;

//***** 全ひずみ変動量の計算 *************************************
//--- ec11の計算 ---
	iii=1; jjj=1;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//xr[i][j]=qrh1[i][j]*zcij(i, j, iii, jjj)+qih1[i][j]*zcij(i, j, iii, jjj);
			//xi[i][j]=qrh1[i][j]*zcij(i, j, iii, jjj)+qih1[i][j]*zcij(i, j, iii, jjj);
			xr[i*ND+j]=qrh1[i*ND+j]*zcij(i, j, iii, jjj, ND)+qih1[i*ND+j]*zcij(i, j, iii, jjj, ND);
			xi[i*ND+j]=qrh1[i*ND+j]*zcij(i, j, iii, jjj, ND)+qih1[i*ND+j]*zcij(i, j, iii, jjj, ND);
		}
	}
	qs=1; rcfft(xrf, xif, xr, xi, ik, s, c, ND, ig);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//ec11[i][j]=xr[i][j];
			ec11[i*ND+j]=xr[i*ND+j];
		}
	}

//--- ec22の計算 ---
	iii=2; jjj=2;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//xr[i][j]=qrh1[i][j]*zcij(i, j, iii, jjj)+qih1[i][j]*zcij(i, j, iii, jjj);
			//xi[i][j]=qrh1[i][j]*zcij(i, j, iii, jjj)+qih1[i][j]*zcij(i, j, iii, jjj);
			xr[i*ND+j]=qrh1[i*ND+j]*zcij(i, j, iii, jjj, ND)+qih1[i*ND+j]*zcij(i, j, iii, jjj, ND);
			xi[i*ND+j]=qrh1[i*ND+j]*zcij(i, j, iii, jjj, ND)+qih1[i*ND+j]*zcij(i, j, iii, jjj, ND);
		}
	}
	qs=1; rcfft(xrf, xif, xr, xi, ik, s, c, ND, ig);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//ec22[i][j]=xr[i][j];
			ec22[i*ND+j]=xr[i*ND+j];
		}
	}

//--- ec12の計算 ---
	iii=1; jjj=2;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//xr[i][j]=qrh1[i][j]*zcij(i, j, iii, jjj)+qih1[i][j]*zcij(i, j, iii, jjj);
			//xi[i][j]=qrh1[i][j]*zcij(i, j, iii, jjj)+qih1[i][j]*zcij(i, j, iii, jjj);
			xr[i*ND+j]=qrh1[i*ND+j]*zcij(i, j, iii, jjj, ND)+qih1[i*ND+j]*zcij(i, j, iii, jjj, ND);
			xi[i*ND+j]=qrh1[i*ND+j]*zcij(i, j, iii, jjj, ND)+qih1[i*ND+j]*zcij(i, j, iii, jjj, ND);
		}
	}
	qs=1; rcfft(xrf, xif, xr, xi, ik, s, c, ND, ig);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//ec12[i][j]=xr[i][j];
			ec12[i*ND+j]=xr[i*ND+j];
		}
	}


//***** 全ひずみ場、応力場、および弾性歪エネルギー場の計算 *************************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//全ひずみ場の計算[式(5.25)]
			//ep11c[i][j]=ec11[i][j]+ep000*c0;
			//ep22c[i][j]=ec22[i][j]+ep000*c0;
			//ep12c[i][j]=ec12[i][j];
			ep11c[i*ND+j]=ec11[i*ND+j]+ep000*c0;
			ep22c[i*ND+j]=ec22[i*ND+j]+ep000*c0;
			ep12c[i*ND+j]=ec12[i*ND+j];

			//弾性応力場の計算[式(5.27)]
			//sig11[i][j]=c11*ec11[i][j]+c12*ec22[i][j]-(c11+2.*c12)*ep000*(ch[i][j]-c0);
			//sig22[i][j]=c12*ec11[i][j]+c11*ec22[i][j]-(c11+2.*c12)*ep000*(ch[i][j]-c0);
			//sig33[i][j]=c12*ec11[i][j]+c12*ec22[i][j]-(c11+2.*c12)*ep000*(ch[i][j]-c0);
			//sig12[i][j]=2.*c44*ec12[i][j];
			sig11[i*ND+j]=c11*ec11[i*ND+j]+c12*ec22[i*ND+j]-(c11+2.0*c12)*ep000*(ch[i*ND+j]-c0);
			sig22[i*ND+j]=c12*ec11[i*ND+j]+c11*ec22[i*ND+j]-(c11+2.0*c12)*ep000*(ch[i*ND+j]-c0);
			sig33[i*ND+j]=c12*ec11[i*ND+j]+c12*ec22[i*ND+j]-(c11+2.0*c12)*ep000*(ch[i*ND+j]-c0);
			sig12[i*ND+j]=2.0*c44*ec12[i*ND+j];

			//弾性ひずみエネルギー場の計算[式(5.28)]
			//Estr1=1.5*(c11+2.0*c12)*ep000*ep000*ch[i][j]*ch[i][j];
			//Estr2=-(c11+2.0*c12)*(ep11c[i][j]+ep22c[i][j])*ep000*ch[i][j];
			//Estr3=0.5*c11*(ep11c[i][j]*ep11c[i][j]+ep22c[i][j]*ep22c[i][j])
            //+c12*ep11c[i][j]*ep22c[i][j]+2.*c44*ep12c[i][j]*ep12c[i][j];
			//Estr[i][j]=Estr1+Estr2+Estr3;
			Estr1=1.5*(c11+2.0*c12)*ep000*ep000*ch[i*ND+j]*ch[i*ND+j];
			Estr2=-(c11+2.0*c12)*(ep11c[i*ND+j]+ep22c[i*ND+j])*ep000*ch[i*ND+j];
			Estr3=0.5*c11*(ep11c[i*ND+j]*ep11c[i*ND+j]+ep22c[i*ND+j]*ep22c[i*ND+j])
                      +c12*ep11c[i*ND+j]*ep22c[i*ND+j]+2.0*c44*ep12c[i*ND+j]*ep12c[i*ND+j];
			Estr[i*ND+j]=Estr1+Estr2+Estr3;
		}
	}


//***** 変位場の計算[式(5.31)] *************************************
//--- u1の計算 ---
	iii=1;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//xr[i][j]=-qrh1[i][j]*zuij(i, j, iii)-qih1[i][j]*zuij(i, j, iii);
			//xi[i][j]=qrh1[i][j]*zuij(i, j, iii)+qih1[i][j]*zuij(i, j, iii);
			xr[i*ND+j]=-qrh1[i*ND+j]*zuij(i, j, iii, ND)-qih1[i*ND+j]*zuij(i, j, iii, ND);
			xi[i*ND+j]= qrh1[i*ND+j]*zuij(i, j, iii, ND)+qih1[i*ND+j]*zuij(i, j, iii, ND);
		}
	}
	qs=1; rcfft(xrf, xif, xr, xi, ik, s, c, ND, ig);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//u1[i][j]=xr[i][j];
			u1[i*ND+j]=xr[i*ND+j];
		}
	}

//--- u2の計算 ---
	iii=2;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//xr[i][j]=-qrh1[i][j]*zuij(i, j, iii)-qih1[i][j]*zuij(i, j, iii);
			//xi[i][j]=qrh1[i][j]*zuij(i, j, iii)+qih1[i][j]*zuij(i, j, iii);
			xr[i*ND+j]=-qrh1[i*ND+j]*zuij(i, j, iii, ND)-qih1[i*ND+j]*zuij(i, j, iii, ND);
			xi[i*ND+j]= qrh1[i*ND+j]*zuij(i, j, iii, ND)+qih1[i*ND+j]*zuij(i, j, iii, ND);
		}
	}
	qs=1; rcfft(xrf, xif, xr, xi, ik, s, c, ND, ig);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//u2[i][j]=xr[i][j];
			u2[i*ND+j]=xr[i*ND+j];
		}
	}

//--- u3の計算 ---
//	iii=3;
//	for(i=0;i<=ndm;i++){
//		for(j=0;j<=ndm;j++){
//			//xr[i][j]=-qrh1[i][j]*zuij(i, j, iii)-qih1[i][j]*zuij(i, j, iii);
//			//xi[i][j]=qrh1[i][j]*zuij(i, j, iii)+qih1[i][j]*zuij(i, j, iii);
//			xr[i*ND+j]=-qrh1[i*ND+j]*zuij(i, j, iii, ND)-qih1[i*ND+j]*zuij(i, j, iii, ND);
//			xi[i*ND+j]= qrh1[i*ND+j]*zuij(i, j, iii, ND)+qih1[i*ND+j]*zuij(i, j, iii, ND);
//		}
//	}
//	qs=1; rcfft(xrf, xif, xr, xi, ik, s, c, ND, ig);
//	for(i=0;i<=ndm;i++){
//		for(j=0;j<=ndm;j++){
//			//u3[i][j]=xr[i][j];
//			u3[i*ND+j]=xr[i*ND+j];
//		}
//	}


//***** 全数値データ保存 ********************************************
	//datsave(ch, Estr, ep11c, ep22c, ep12c,
	//				sig11, sig22, sig33, sig12, u1, u2, ND);	//save
	datsave_paraview(ch, Estr, ep11c, ep22c, ep12c,
					sig11, sig22, sig33, sig12, u1, u2, ND);	//save
	//if(keypress()){return 0;}//キー待ち状態
	return 0;
}


//******* Sin, Cos のテーブルおよびビット反転テーブルの設定 ***************
void table(double *s, double *c, int *ik, int ND, int ig)
{
	int i, j, mc, mn;
	double q;
	int nd=ND, nd2=ND/2;

	q=2.0*PI/nd;
	for(i=0;i<=nd2-1;i++){ c[i]=cos(q*i); s[i]=sin(q*i); }//Sin, Cos のテーブル

	ik[0]=0; mn=nd2; mc=1;
	for(i=1;i<=ig;i++){
		for(j=0;j<=mc-1;j++){ ik[j+mc]=ik[j]+mn; }	//ビット反転テーブル
		mc*=2; mn/=2; 
	}
}

//********** １次元高速フーリエ変換 **************************************
void fft(double *xrf, double *xif, double *s, double *c, int ND, int ig)
{
	int lg, lf, mf, nf, n2, ka, kb, ix;
	double tj, tr;
	int nd2=ND/2;

	lg=1;
	for(lf=1;lf<=ig;lf++){
		n2=nd2/lg;
		for(mf=1;mf<=lg;mf++){
			for(nf=0;nf<=n2-1;nf++){
				ix=nf*lg;
				ka=nf+2*n2*(mf-1);
				kb=ka+n2;
				tr=xrf[ka]-xrf[kb];
				tj=xif[ka]-xif[kb];
				xrf[ka]+=xrf[kb];
				xif[ka]+=xif[kb];
				xrf[kb]=tr*c[ix]-tj*float(qs)*s[ix];
				xif[kb]=tj*c[ix]+tr*float(qs)*s[ix];
			}
		}
		lg*=2;
	}//lf
}

//************ ２次元高速フーリエ変換 ***********************************
void rcfft(double *xrf, double *xif, double *xr, double *xi, int *ik, double *s, double *c, int ND, int ig)
{
	int i, ic, ir, j;
	int nd=ND, ndm=ND-1;

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

	if(qs > 0){return;}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//xr[i][j]=xr[i][j]/nd/nd;  xi[i][j]=xi[i][j]/nd/nd;
			xr[i*ND+j]=xr[i*ND+j]/nd/nd;
			xi[i*ND+j]=xi[i*ND+j]/nd/nd;
		}
	}

}

//*** Zcij [式(5.26)の計算] ****************************************
double zcij(int i0, int j0, int iii, int jjj, int ND)
{
	//int i, j, k, m, n, p, q;
	int m, n;
	int ii, jj;
	double nx, ny, nz, alnn;
	double nec[4];
	double a11, a22, a33, a12, a13, a23, a21, a31, a32;
	double b11, b22, b33, b12, b13, b23, b21, b31, b32;
	double det1, zij;
	double om[4][4];
	int nd=ND, nd2=ND/2;

	if(i0<=nd2-1) {ii=i0;}
	if(i0>=nd2) {ii=i0-nd;}
	if(j0<=nd2-1) {jj=j0;}
	if(j0>=nd2) {jj=j0-nd;}
	alnn=sqrt((double)ii*(double)ii+(double)jj*(double)jj);
	if(alnn==0.0){alnn=1.0;}
	nec[1]=nx=(double)ii/alnn;
	nec[2]=ny=(double)jj/alnn;
	nec[3]=nz=0.0;

  //for(i=1;i<=3;i++){ for(j=1;j<=3;j++){ om[i][j]=0.; } }

	a11=cec[1][1][1][1]*nx*nx+cec[1][2][1][2]*ny*ny+cec[1][3][1][3]*nz*nz;
	a22=cec[1][2][1][2]*nx*nx+cec[2][2][2][2]*ny*ny+cec[2][3][2][3]*nz*nz;
	a33=cec[3][1][3][1]*nx*nx+cec[2][3][2][3]*ny*ny+cec[3][3][3][3]*nz*nz;
	a12=(cec[1][1][2][2]+cec[1][2][1][2])*nx*ny;
	a23=(cec[2][2][3][3]+cec[2][3][2][3])*ny*nz;
	a31=(cec[3][3][1][1]+cec[3][1][3][1])*nx*nz;
	a21=a12;
	a32=a23;
	a13=a31;

	b11=a22*a33-a23*a32;
	b22=a11*a33-a13*a31;
	b33=a11*a22-a12*a21;
	b12=-(a21*a33-a23*a31);
	b23=-(a11*a32-a12*a31);
	b31=-(a22*a31-a21*a32);
	b21=b12;
	b32=b23;
	b13=b31;

	det1=a11*a22*a33+a12*a23*a31+a13*a32*a21-a13*a31*a22-a11*a23*a32-a33*a12*a21;
	if(det1==0.0){det1=1.0;}

	om[1][1]=b11/det1;
	om[2][2]=b22/det1;
	om[3][3]=b33/det1;
	om[1][2]=om[2][1]=b12/det1;
	om[2][3]=om[3][2]=b23/det1;
	om[3][1]=om[1][3]=b31/det1;

	zij=0.0;
	for(m=1;m<=3;m++) {
		for(n=1;n<=3;n++) {
			zij=zij+0.5*(sigma[m][n]*nec[jjj]*nec[n]*om[m][iii]
						+sigma[m][n]*nec[iii]*nec[n]*om[m][jjj]);
		}
   }
	return(zij);
}

//*** Zuij  [式(5.30)の計算] ****************************************
double zuij(int i0, int j0, int iii, int ND)
{
	//int i, j, k, m, n, p, q;
	int m, n;
	int ii, jj;
	double nx, ny, nz, alnn;
	double nec[4];
	double a11, a22, a33, a12, a13, a23, a21, a31, a32;
	double b11, b22, b33, b12, b13, b23, b21, b31, b32;
	double det1, zij;
	double om[4][4];
	int nd=ND, nd2=ND/2;

	if(i0<=nd2-1) {ii=i0;}
	if(i0>=nd2) {ii=i0-nd;}
	if(j0<=nd2-1) {jj=j0;}
	if(j0>=nd2) {jj=j0-nd;}
	alnn=sqrt((double)ii*(double)ii+(double)jj*(double)jj);
	if(alnn==0.0){alnn=1.0;}
	nec[1]=nx=(double)ii/alnn;
	nec[2]=ny=(double)jj/alnn;
	nec[3]=nz=0.0;

  //for(i=1;i<=3;i++){ for(j=1;j<=3;j++){ om[i][j]=0.; } }

	a11=cec[1][1][1][1]*nx*nx+cec[1][2][1][2]*ny*ny+cec[1][3][1][3]*nz*nz;
	a22=cec[1][2][1][2]*nx*nx+cec[2][2][2][2]*ny*ny+cec[2][3][2][3]*nz*nz;
	a33=cec[3][1][3][1]*nx*nx+cec[2][3][2][3]*ny*ny+cec[3][3][3][3]*nz*nz;
	a12=(cec[1][1][2][2]+cec[1][2][1][2])*nx*ny;
	a23=(cec[2][2][3][3]+cec[2][3][2][3])*ny*nz;
	a31=(cec[3][3][1][1]+cec[3][1][3][1])*nx*nz;
	a21=a12;
	a32=a23;
	a13=a31;

	b11=a22*a33-a23*a32;
	b22=a11*a33-a13*a31;
	b33=a11*a22-a12*a21;
	b12=-(a21*a33-a23*a31);
	b23=-(a11*a32-a12*a31);
	b31=-(a22*a31-a21*a32);
	b21=b12;
	b32=b23;
	b13=b31;

	det1=a11*a22*a33+a12*a23*a31+a13*a32*a21-a13*a31*a22-a11*a23*a32-a33*a12*a21;
	if(det1==0.0){det1=1.0;}

	om[1][1]=b11/det1;
	om[2][2]=b22/det1;
	om[3][3]=b33/det1;
	om[1][2]=om[2][1]=b12/det1;
	om[2][3]=om[3][2]=b23/det1;
	om[3][1]=om[1][3]=b31/det1;

	zij=0.0;
	for(m=1;m<=3;m++) {
		for(n=1;n<=3;n++) {
			zij=zij+sigma[m][n]*nec[n]*om[m][iii];
		}
	}
	zij=-zij/alnn;
	return(zij);
}

//************ 各種の場のデータの保存 *******************************
void datsave(double *ch, double *Estr, 
			 double *ep11c, double *ep22c, double *ep12c,
			 double *sig11, double *sig22, double *sig33, double *sig12, 
			 double *u1, double *u2, int ND)
{
	FILE		*stream;//ストリームのポインタ設定
	int 		i, j;//整数
	int ndm=ND-1;

	stream = fopen("el_field.dat", "w");//書き込み先のファイルを上書き方式でオープン

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			fprintf(stream, "%e  ", ch[i*ND+j]);//各種場のデータ保存
			fprintf(stream, "%e  ", Estr[i*ND+j]);
			fprintf(stream, "%e  ", ep11c[i*ND+j]);
			fprintf(stream, "%e  ", ep22c[i*ND+j]);
			fprintf(stream, "%e  ", ep12c[i*ND+j]);
			fprintf(stream, "%e  ", sig11[i*ND+j]);
			fprintf(stream, "%e  ", sig22[i*ND+j]);
			fprintf(stream, "%e  ", sig33[i*ND+j]);
			fprintf(stream, "%e  ", sig12[i*ND+j]);
			fprintf(stream, "%e  ", u1[i*ND+j]);
			fprintf(stream, "%e  ", u2[i*ND+j]);
		}
	}
	fprintf(stream, "\n");//改行の書き込み
	fclose(stream);//ファイルをクローズ
}

void datsave_paraview(double *ch, double *Estr, 
			 double *ep11c, double *ep22c, double *ep12c,
			 double *sig11, double *sig22, double *sig33, double *sig12, 
			 double *u1, double *u2, int ND)
{
	FILE	*fp;
	char	fName[256];
	int 	i, j;
	int ndm=ND-1;
	
	iout = iout + 1;
	printf("pf_result%06d.vtk \n",iout);
	sprintf(fName,"pf_result%06d.vtk",iout);
	fp = fopen(fName, "w");
	fprintf(fp,"# vtk DataFile Version 3.0 \n");
	fprintf(fp,"output.vtk \n");
	fprintf(fp,"ASCII \n");
	fprintf(fp,"DATASET STRUCTURED_POINTS \n");
	fprintf(fp,"DIMENSIONS %5d %5d %5d \n",(ndm+1),(ndm+1),1);
	fprintf(fp,"ORIGIN 0.0 0.0 0.0 \n");
	fprintf(fp,"ASPECT_RATIO 1 1 1 \n");
	fprintf(fp,"POINT_DATA %11d \n",((ndm+1)*(ndm+1)*1));
	fprintf(fp,"SCALARS concentration float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			fprintf(fp,"%10.6f\n", ch[i*ND+j]);//各種場のデータ保存
		}
	}
	fprintf(fp,"SCALARS Elastic_strain_energy_density float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			fprintf(fp,"%10.6f\n", Estr[i*ND+j]);
		}
	}
	fprintf(fp,"SCALARS Strain_c11 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			fprintf(fp,"%10.6f\n", ep11c[i*ND+j]);
		}
	}
	fprintf(fp,"SCALARS Strain_c22 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			fprintf(fp,"%10.6f\n", ep22c[i*ND+j]);
		}
	}
	fprintf(fp,"SCALARS Strain_c12 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			fprintf(fp,"%10.6f\n", ep12c[i*ND+j]);
		}
	}
	fprintf(fp,"SCALARS Elastic_stress_11 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			fprintf(fp,"%10.6f\n", sig11[i*ND+j]);
		}
	}
	fprintf(fp,"SCALARS Elastic_stress_22 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			fprintf(fp,"%10.6f\n", sig22[i*ND+j]);
		}
	}
	fprintf(fp,"SCALARS Elastic_stress_33 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			fprintf(fp,"%10.6f\n", sig33[i*ND+j]);
		}
	}
	fprintf(fp,"SCALARS Elastic_stress_12 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			fprintf(fp,"%10.6f\n", sig12[i*ND+j]);
		}
	}
	fprintf(fp,"SCALARS Displacement_u1 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			fprintf(fp,"%10.6f\n", u1[i*ND+j]);
		}
	}
	fprintf(fp,"SCALARS Displacement_u2 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			fprintf(fp,"%10.6f\n", u2[i*ND+j]);
		}
	}
	fprintf(fp,"VECTORS vectors float \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			fprintf(fp,"%10.6f %10.6f %10.6f\n", u1[i*ND+j], u2[i*ND+j], 0.0); // 変位場を矢印にて描画
		}
	}
	fclose(fp);
}
//************ 場のデータの読み込み *****************************************
void datin(double *ch, int ND)
{
	FILE		*datin0;//ストリームのポインタ設定
	int 		i, j;//整数
	double c00;//場の平均値
	int nd=ND, ndm=ND-1;

	datin0 = fopen("ph.dat", "r");//読込み元のファイルをオープン
	c00=0.0;//場の平均値の初期値
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			fscanf(datin0, "%lf", &ch[i*ND+j]);//場のデータ読込
			c00+=ch[i*ND+j];
		}
	}
	c0=c00/nd/nd;//場の平均値
	printf("c0=  %f  \n", c0);

	fclose(datin0);//ファイルをクローズ
}
