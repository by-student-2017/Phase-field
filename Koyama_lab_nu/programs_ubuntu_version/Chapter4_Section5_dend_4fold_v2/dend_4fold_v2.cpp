#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

#define DRND(x) ((double)(x)/RAND_MAX*rand())//乱数の関数設定

//#define NDP 351							//差分計算における計算領域一辺の分割数＋１

	//int nd=NDP-1, ndm=NDP-2, nd2=(NDP-1)/2;//計算領域の一辺の差分ブロック数、nd-1を定義、nd/2を定義
	double delt;						//時間刻み
	double PI=3.14159;					//円周率
	double RR=8.3145;					//ガス定数
	double Tini, Tm, time1;				//初期温度、融点、計算カウント数
	//double s1h[NDP][NDP], Th[NDP][NDP];	//フェーズフィールド、温度場
	int Nstep, iout;

	void ini000(double *s1h, double *Th, int NDP);	//初期場の設定サブル－チン
	void datsave(double *s1h, double *Th, int NDP);	//デ－タ保存サブル－チン
	void datsave_paraview(double *s1h, double *Th, int NDP);//デ－タ保存サブル－チン
	

//******* メインプログラム ******************************************
int main(void)
{
	int NDP;						//差分計算における計算領域一辺の分割数＋１
	int nd, ndm, nd2;				//計算領域の一辺の差分ブロック数、nd-1を定義、nd/2を定義
	
	//double s1h2[NDP][NDP], Th2[NDP][NDP];	//場の補助配列
	double s1kai, s1kais;			//ポテンシャル
	int   i, j, k, l;				//整数
	int   ip, im, jp, jm;			//整数
	double time1max;				//計算カウント数の最大値（計算終了カウント）

	double s1;						//フェーズフィールド
	double s1ip, s1im, s1jp, s1jm;	//s1の左右上下のフェーズフィールド
	double s1ipjp, s1ipjm, s1imjp, s1imjm;//s1の斜め横のフェーズフィールド
	double s1ddtt;					//s1の時間変化量（発展方程式の左辺）

	double TT;						//温度場
	double Tip, Tim, Tjp, Tjm;		//TTの左右上下の温度
	double Tddtt;					//温度の時間変化量（熱拡散方程式の左辺）

	double ep, ep1p, ep2p;			//勾配エネルギー係数関連の変数
	double dx_s1, dy_s1, l_dxdy;	//フェーズフィールドの空間微分関連の変数
	double dxx_s1, dyy_s1, dxy_s1;	//フェーズフィールドの空間微分関連の変数

	double al;				//計算領域の１辺の長さ
	double dx, dy;			//差分格子サイズ(x方向)
	double gamma;			//界面エネルギー密度
	double delta;			//界面幅
	double ram;				//λ
	double bbb;				//界面幅に関係する係数
	double j_fold;			//異方性モード
	double astre;			//異方性強度
	double cndct;			//熱伝導率
	double speht;			//比熱
	double rlate;			//潜熱
	double skine;			//界面カイネティック係数
	double th0;				//優先成長方向の角度
	double th;				//界面の法線方向の角度
	double aaa, aaac;		//勾配エネルギー係数0
	double www, wwwc;		//ペナルティー項のエネルギー障壁
	double pmobi, pmobic;	//フェーズフィールドのモビリティ
	double dtp;				//時間刻み
	double dtt;				//時間刻み
	double anois;			//ノイズの振幅
	double dF;				//駆動力
	double dami1, dami2;	//場の計算をスキップする条件設定に使用する変数

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
	NDP     = int(data[0]);
	delt    = data[1];
	dx      = data[2];
	dy      = data[3];
	gamma   = data[4];
	ram     = data[5];
	j_fold  = data[6];
	astre   = data[7];
	cndct   = data[8];
	speht   = data[9];
	rlate   = data[10];
	Tm      = data[11];
	Tini    = data[12];
	skine   = data[13];
	th0     = data[14];
	aaac    = data[15];
	wwwc    = data[16];
	pmobic  = data[17];
	anois   = data[18];
	time1max= int(data[19]);
	Nstep   = int(data[20]);
	printf("---------------------------------\n");
	//
	nd=NDP-1;
	ndm=NDP-2;
	nd2=(NDP-1)/2;
	//
	double *s1h  = (double *)malloc(sizeof(double)*( NDP*NDP + NDP ));	//フェーズフィールド
	double *Th   = (double *)malloc(sizeof(double)*( NDP*NDP + NDP ));	//温度場
	double *s1h2 = (double *)malloc(sizeof(double)*( NDP*NDP + NDP ));
	double *Th2  = (double *)malloc(sizeof(double)*( NDP*NDP + NDP ));
	//
	//printf("dtemp(1.0)=  ");	scanf(" %lf",&dtemp);	//dtemp=1.0;
	//printf("DELT(1.5)=  "); scanf(" %lf",&delt);	//delt=1.5;

	//dx=dy=30.0e-9;				//差分格子サイズ(x方向)[m]
	delta=3.0*dx;      		//界面幅[m]
	//dx=dy=20.0e-9;			//差分格子サイズ(x方向)[m]
	//delta=4.0*dx;    		//界面幅[m]
	al=dx*(double)nd;			//計算領域[m]
	//gamma=0.37;       		//界面エネルギー[J/m^2]
	//ram=0.1;        			//λ
	bbb=2.0*log((1.0+(1.0-2.0*ram))/(1.0-(1.0-2.0*ram)))/2.0;//界面幅に関係する係数
	//j_fold=4.0;         	//異方性モード

	//astre=0.005;       	//異方性強度
	//astre=0.01;       	//異方性強度
	//astre=0.03;       		//異方性強度(図4.8)
	//astre=0.05;       	//異方性強度

//    << pure Ni >>
	//cndct=84.01;      		//熱伝導率[W/mK]
	//speht=5.42e+06;   		//比熱[J/Km^3]
	//rlate=2.350e+09;  		//潜熱[J/m^3]
	//Tm=1728.0;      		//融点[K]
	//Tini=1511.2;     		//初期温度[K]
	//Tini=1424.5;     		//初期温度[K]
	//skine=2.0;        		//界面カイネティック係数[m/Ks]
	//th0=0.0;         		//優先成長方向の角度

//---- 4.5.2項参照 ----------------------------------------------------
	//aaa=sqrt(3.0*delta*gamma/bbb);         	//勾配エネルギー係数0
	aaa=sqrt(aaac*delta*gamma/bbb);         	//勾配エネルギー係数
	//www=6.0*gamma*bbb/delta;               	//ペナルティー項のエネルギー障壁
	www=wwwc*gamma*bbb/delta;               	//ペナルティー項のエネルギー障壁
	//pmobi=bbb*Tm*skine/(3.0*delta*rlate); 	//フェーズフィールドのモビリティ
	pmobi=bbb*Tm*skine/(pmobic*delta*rlate); 	//フェーズフィールドのモビリティ

	dtp=dx*dx/(5.0*pmobi*aaa*aaa);	//時間きざみ
	dtt=dx*dx/(5.0*cndct/speht);	//時間きざみ
	if(dtp>dtt){delt=dtt;} else{delt=dtp;}
	printf("delt= %e \n", delt);
//-----------------------------------------------------------------

	//anois=0.1;	//ノイズの振幅

	time1=0.0;		//計算時間の初期値
	//time1max=1.0+1.0e+08;	//計算時間の最大値

//*** 初期濃度場の設定と描画Window表示 *****************************************

	ini000(s1h, Th, NDP);//初期場の設定

//**** シミュレーションスタート ******************************
//Nstep = 100;
iout = -1;
start: ;

	if((((int)(time1) % Nstep)==0)){datsave(s1h, Th, NDP);}		//一定繰返しカウント毎に場を保存
	if((((int)(time1) % Nstep)==0)){datsave_paraview(s1h, Th, NDP);}//一定繰返しカウント毎に場を保存

//****** フェーズフィールドおよび温度場の時間発展  **************
	for(i=0;i<=nd;i++){
		for(j=0;j<=nd;j++){
			ip=i+1; im=i-1; jp=j+1; jm=j-1;
			if(i==nd){ip=ndm;}  if(i==0){im=1;}
			if(j==nd){jp=ndm;}  if(j==0){jm=1;}

			//s1=s1h[i][j];//フェーズフィールド
			//s1ip=s1h[ip][j];  s1im=s1h[im][j];  s1jp=s1h[i][jp];  s1jm=s1h[i][jm];
			//s1ipjp=s1h[ip][jp]; s1ipjm=s1h[ip][jm];  s1imjp=s1h[im][jp];  s1imjm=s1h[im][jm];
			s1=s1h[i*NDP+j];//フェーズフィールド
			s1ip=s1h[ip*NDP+j];  s1im=s1h[im*NDP+j];  s1jp=s1h[i*NDP+jp];  s1jm=s1h[i*NDP+jm];
			s1ipjp=s1h[ip*NDP+jp]; s1ipjm=s1h[ip*NDP+jm];  s1imjp=s1h[im*NDP+jp];  s1imjm=s1h[im*NDP+jm];

			//TT=Th[i][j];//温度場
			//Tip=Th[ip][j];  Tim=Th[im][j];  Tjp=Th[i][jp];  Tjm=Th[i][jm];
			TT=Th[i*NDP+j];//温度場
			Tip=Th[ip*NDP+j];  Tim=Th[im*NDP+j];  Tjp=Th[i*NDP+jp];  Tjm=Th[i*NDP+jm];

//----- 計算をスキップする場合の判断 ----------------------------------------------
			dami1=fabs(s1+s1ip+s1im+s1jp+s1jm);  dami2=fabs(TT+Tip+Tim+Tjp+Tjm-5.0*Tini);
			if( (dami1<=1.0e-20)&&(dami2<=1.0e-20) ){
				//s1h2[i][j]=s1h[i][j];  Th2[i][j]=Th[i][j];
				s1h2[i*NDP+j]=s1h[i*NDP+j];  Th2[i*NDP+j]=Th[i*NDP+j];
				goto dami;
			}
//---------------------------------------------------------------------------------

			dx_s1=(s1ip-s1im)/2.0/dx;  				//フェーズフィールドの空間１階微分
			dy_s1=(s1jp-s1jm)/2.0/dy;
			dxx_s1=(s1ip+s1im-2.0*s1)/dx/dx;		//フェーズフィールドの空間２階微分
			dyy_s1=(s1jp+s1jm-2.0*s1)/dy/dy;
			dxy_s1=(s1ipjp+s1imjm-s1imjp-s1ipjm)/4.0/dx/dy;
			th=atan(dy_s1/(dx_s1+1.0e-20));			//界面の法線方向の角度[式(4.24)]

			ep=aaa*(1.0+astre*cos(j_fold*(th-th0)));	//勾配エネルギー係数の平方根[式(4.23)]
			ep1p=-aaa*astre*j_fold*sin(j_fold*(th-th0));	//epの角度による１階微分
			ep2p=-aaa*astre*j_fold*j_fold*cos(j_fold*(th-th0));	//epの角度による２階微分

			s1kais=-ep*ep*(dxx_s1+dyy_s1)
					-ep*ep1p*((dyy_s1-dxx_s1)*sin(2.0*th)+2.0*dxy_s1*cos(2.0*th))
					 +0.5*(ep1p*ep1p+ep*ep2p)*(2.0*dxy_s1*sin(2.0*th)
					-dxx_s1-dyy_s1-(dyy_s1-dxx_s1)*cos(2.0*th));
			//勾配ポテンシャル[式(4.29)の一部]

			dF=15.0/(2.0*www)*rlate*(TT-Tm)/Tm*s1*(1.0-s1);//化学的駆動力[式(4.29)の一部]
			s1kai=4.0*www*s1*(1.0-s1)*(0.5-s1+dF+anois*(DRND(1)-0.5));
			//化学ポテンシャル[式(4.29)の一部]

			s1ddtt=-pmobi*(s1kai+s1kais);	//フェーズフィールドの発展方程式[式(4.25)]

			//s1h2[i][j]=s1+s1ddtt*delt;	//フェーズフィールドの時間発展(陽解法)
			//s1h2[i][j]=s1+s1ddtt*delt+anois*(DRND(1)-0.5)*s1*(1.0-s1);
			//if(s1h2[i][j]>=1.0){s1h2[i][j]=1.0;}//フェーズフィールドの変域補正
			//if(s1h2[i][j]<=0.0){s1h2[i][j]=0.0;}
			s1h2[i*NDP+j]=s1+s1ddtt*delt;	//フェーズフィールドの時間発展(陽解法)
			if(s1h2[i*NDP+j]>=1.0){s1h2[i*NDP+j]=1.0;}//フェーズフィールドの変域補正
			if(s1h2[i*NDP+j]<=0.0){s1h2[i*NDP+j]=0.0;}

			Tddtt=( cndct*( (Tip+Tim-2.0*TT)/dx/dx+(Tjp+Tjm-2.0*TT)/dy/dy )
      	                    +30.0*s1*(1.0-s1)*s1*(1.0-s1)*rlate*s1ddtt )/speht;
			//熱拡散方程式[式(4.30)]
			//Th2[i][j]=Th[i][j]+Tddtt*delt;		//温度場の時間発展(陽解法)
			Th2[i*NDP+j]=Th[i*NDP+j]+Tddtt*delt;		//温度場の時間発展(陽解法)

			dami:;
		}
	}

//**** 場の補助配列を主配列に移動 ************************************
	for(i=0;i<=nd;i++){
	    //for(j=0;j<=nd;j++){  s1h[i][j]=s1h2[i][j]; Th[i][j]=Th2[i][j]; }
		for(j=0;j<=nd;j++){  s1h[i*NDP+j]=s1h2[i*NDP+j]; Th[i*NDP+j]=Th2[i*NDP+j]; }
	}
//*********************************************************************

	//if(keypress()){return 0;}//キー待ち状態
	time1=time1+1.;  if(time1<time1max){goto start;}//時間増加および最大時間に到達したかどうかの判断

	end:;
  return 0;

}

//************ 初期場の設定サブルーチン *************
void ini000(double *s1h, double *Th, int NDP)
{
	int i, j;
 	srand(time(NULL)); // 乱数初期化
	int nd=NDP-1, ndm=NDP-2, nd2=(NDP-1)/2;

	for(i=0;i<=nd;i++){
		for(j=0;j<=nd;j++){
			//s1h[i][j]=0.0;
			//if((i<10.)&&(j<10.)){s1h[i][j]=0.9;}
			//if((i*i+j*j)<20.){s1h[i][j]=0.9;}
			s1h[i*NDP+j]=0.0;
			if((i*i+j*j)<20.0){s1h[i*NDP+j]=0.9;}
		}
	}

	for(i=0;i<=nd;i++){
		for(j=0;j<=nd;j++){
			//Th[i][j]=Tini+s1h[i][j]*(Tm-Tini);
			Th[i*NDP+j]=Tini+s1h[i*NDP+j]*(Tm-Tini);
		}
	}

}

//************ データ保存サブルーチン *******************************
void datsave(double *s1h, double *Th, int NDP)
{
	FILE		*stream;	//ストリームのポインタ設定
	int 		i, j;			//整数
	int nd=NDP-1, ndm=NDP-2, nd2=(NDP-1)/2;

	stream = fopen("test.dat", "a");	//書き込む先のファイルを追記方式でオープン
	fprintf(stream, "%f\n", time1);		//繰返しカウント数の保存
	for(i=0;i<=nd;i++){
		for(j=0;j<=nd;j++){
			//fprintf(stream, "%e  %e  ", s1h[i][j], Th[i][j]);//局所場の保存
			fprintf(stream, "%e  %e  ", s1h[i*NDP+j], Th[i*NDP+j]);//局所場の保存
		}
	}
	fprintf(stream, "\n");	//改行の書き込み
	fclose(stream);					//ファイルをクローズ
}

void datsave_paraview(double *s1h, double *Th, int NDP)
{
	FILE	*fp;
	char	fName[256];
	int 	i, j;
	int nd=NDP-1, ndm=NDP-2, nd2=(NDP-1)/2;
	
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
	fprintf(fp,"SCALARS Phase_field float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			//fprintf(stream, "%e  %e  ", s1h[i][j], Th[i][j]);//局所場の保存
			//fprintf(fp,"%10.6f\n", s1h[i][j]);
			fprintf(fp,"%10.6f\n", s1h[i*NDP+j]);
		}
	}
	fprintf(fp,"SCALARS Temperature float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			//fprintf(stream, "%e  %e  ", s1h[i][j], Th[i][j]);//局所場の保存
			//fprintf(fp,"%10.6f\n", Th[i][j]);
			fprintf(fp,"%10.6f\n", Th[i*NDP+j]);
		}
	}
	fclose(fp);
}
