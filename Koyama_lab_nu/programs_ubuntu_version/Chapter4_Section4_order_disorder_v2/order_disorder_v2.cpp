#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

#define DRND(x) ((double)(x)/RAND_MAX*rand())//乱数の関数設定

//#define ND 400						//差分計算における計算領域一辺の分割数

	//int nd=ND, ndm=ND-1; 			//計算領域の一辺の差分分割数(差分ブロック数)、ND-1を定義
	//int nd2=ND/2;				 	//ND/2を定義使用
	double PI=3.14159;				//円周率
	double rr=8.3145;				//ガス定数
	double time1;					//時間
	int iout;

	//double s1h[ND][ND];				//長範囲規則度場

	void ini000(double *s1h, int ND);	//初期場の設定サブル－チン
	void datsave(double *s1h, int ND);	//デ－タ保存サブル－チン
	void datsave_paraview(double *s1h, int ND);//デ－タ保存サブル－チン

//******* メインプログラム ******************************************
int main(void)
{
	int ND, nd, ndm, nd2;
	
	double s1;						//長範囲規則度場
	//double s1h2[ND][ND];			//場の補助配列
	double s1k_chem, s1k_su;		//ポテンシャル
	double sum11;					//s1の空間積分
	double s1ddtt;					//s1の時間変化量（発展方程式の左辺）

	int   i, j, k, l, ii, jj, kk, iii, jjj;	//整数
	int   p, q, m, n;						//整数
	int   ip, im, jp, jm, Nstep;			//整数
	double al, temp, delt;					//計算領域、時間、時間きざみ
	double time1max;						//最大時間（計算を止める際に使用）
	double b1, vm0, atom_n;					//規格化長さ、モル体積、単位胞内の原子数
	double smob;							//結晶変態の緩和係数

	double AA0, AA1, AA2, AA3;				//ギズブエネルギー内の係数
	double AA0e;
	double a1_c, b1_c, c1_c;				//格子定数
	double kappa_s1;						//勾配エネルギ－係数
	double kappa_s1c;
	//double ds_fac;							//結晶変態の揺らぎ係数

//****** 計算条件および物質定数の設定 ****************************************
	printf("---------------------------------\n");
	printf("read parameters form parameters.txt\n");
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
	temp    = data[2];
	al      = data[3];	// [nm]
	smob    = data[4];
	AA0e    = data[5];
	kappa_s1c=data[6];
	vm0     = data[7];
	time1max= int(data[8]);
	Nstep   = int(data[9]);
	printf("---------------------------------\n");
	
	//
	nd=ND, ndm=ND-1; 			//計算領域の一辺の差分分割数(差分ブロック数)、ND-1を定義
	nd2=ND/2;				 	//ND/2を定義使用
	//
	double *s1h  = (double *)malloc(sizeof(double)*( ND*ND + ND ));	//長範囲規則度場
	double *s1h2 = (double *)malloc(sizeof(double)*( ND*ND + ND ));	//場の補助配列
	//
	//printf("DELT(1.0)=  ");	scanf(" %lf",&delt);//時間きざみ入力	//delt=1.0;

	//temp=500.0;					//温度(K)
	//al=1000.0*1.0E-09;			//計算領域[m]
	al=al*1.0E-09;				//計算領域[m]
	b1=al/nd;					//差分ブロックの長さ

	time1=0.0;					//初期計算カウント数の設定
	//time1max=1.0+1.0e+07;		//最大計算カウント数の設定

	//smob=1.0;					//モビリティー（変態の緩和係数）

	//AA0=200.0/rr/temp;			//規則-不規則変態の化学的駆動力
	AA0=AA0e/rr/temp;			//規則-不規則変態の化学的駆動力

	//kappa_s1=5.0e-15/rr/temp/b1/b1;//勾配エネルギ－係数
	kappa_s1=kappa_s1c/rr/temp/b1/b1;//勾配エネルギ－係数

	//a1_c=b1_c=c1_c=3.563E-10;	//格子定数
	//atom_n=4.0;  vm0=6.02E23*a1_c*b1_c*c1_c/atom_n;	//モル体積の計算（fccを仮定）

//*** 初期場の設定 ***************

	ini000(s1h, ND);

//**** シミュレーションスタート ******************************
//Nstep = 10;
iout = -1;
start: ;

	//if(time1<=200.){Nstep=10;} else{Nstep=100;}		//データ保存する時間間隔の変更
	if((((int)(time1) % Nstep)==0)) {datsave(s1h, ND);} 	//一定繰返しカウント毎に組織データを保存
	if((((int)(time1) % Nstep)==0)) {datsave_paraview(s1h, ND);} 	//一定繰返しカウント毎に組織データを保存

//******  ポテンシャルの計算 ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ip=i+1; im=i-1; jp=j+1; jm=j-1;
			if(i==ndm){ip=0;}  if(i==0){im=ndm;}
			if(j==ndm){jp=0;}  if(j==0){jm=ndm;}

			//s1=s1h[i][j];
			//s1k_su=-kappa_s1*(s1h[ip][j]+s1h[im][j]+s1h[i][jp]+s1h[i][jm]-4.0*s1);
			s1=s1h[i*ND+j];
			s1k_su=-kappa_s1*(s1h[ip*ND+j]+s1h[im*ND+j]+s1h[i*ND+jp]+s1h[i*ND+jm]-4.0*s1);
															//勾配ポテンシャル[式(4.16)]
			s1k_chem=-4.0*AA0*s1*(1.0-s1*s1);				//化学ポテンシャルの計算[式(4.14)]

			s1ddtt=-smob*(s1k_chem+s1k_su);					//場の時間発展の計算[式(4.17)]
			//s1h2[i][j]=s1h[i][j]+s1ddtt*delt;				//陽解法
			s1h2[i*ND+j]=s1h[i*ND+j]+s1ddtt*delt;				//陽解法

			//if(s1h2[i][j]>=1.0){s1h2[i][j]=1.0;}  
			//if(s1h2[i][j]<=-1.0){s1h2[i][j]=-1.0;}			//sの変域(-1<=s<=1)の補正
			if(s1h2[i*ND+j]>= 1.0){s1h2[i*ND+j]= 1.0;}  
			if(s1h2[i*ND+j]<=-1.0){s1h2[i*ND+j]=-1.0;}			//sの変域(-1<=s<=1)の補正
		}
	}

	//for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){  s1h[i][j]=s1h2[i][j]; } }
	for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){  s1h[i*ND+j]=s1h2[i*ND+j]; } }
	//補助配列から主配列にデータコピー

	//if(keypress()){return 0;}				//キー待ち状態

	time1=time1+1.0;								//計算カウント数の加算
	if(time1<time1max){goto start;}	//最大カウント数に到達したかどうかの判断

end:;
  return 0;
}

//************ 初期波設定サブル－チン *******************************
void ini000(double *s1h, int ND)
{
	int i, j;
	srand(time(NULL)); // 乱数初期化
	int ndm=ND-1;

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//s1h[i][j]=0.01*(2.0*DRND(1.0)-1.0);//場を最大±1%の乱数にて設定
			s1h[i*ND+j]=0.01*(2.0*DRND(1.0)-1.0);//場を最大±1%の乱数にて設定
		}
	}
}

//************ データの保存サブルーチン *******************************
void datsave(double *s1h, int ND)
{
	FILE		*stream;	//ストリームのポインタ設定
	int 		i, j;			//整数
	int ndm=ND-1;

	stream = fopen("test.dat", "a");	//書き込み先のファイルを追記方式でオープン

	fprintf(stream, "%e\n", time1);		//繰返しカウント数の保存
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//fprintf(stream, "%e  ", s1h[i][j]);//場のデータ保存
			fprintf(stream, "%e  ", s1h[i*ND+j]);//場のデータ保存
		}
	}
	fprintf(stream, "\n");	//改行の書き込み
	fclose(stream);					//ファイルをクローズ
}

void datsave_paraview(double *s1h, int ND)
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
	fprintf(fp,"DIMENSIONS %5d %5d %5d \n",(ndm+1),(ndm+1),(1));
	fprintf(fp,"ORIGIN 0.0 0.0 0.0 \n");
	fprintf(fp,"ASPECT_RATIO 1 1 1 \n");
	fprintf(fp,"POINT_DATA %16d \n",((ndm+1)*(ndm+1)*1));
	fprintf(fp,"SCALARS Phase_field float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			//fprintf(stream, "%e  ", s1h[i][j]);//場のデータ保存
			//fprintf(fp,"%10.6f\n", s1h[i][j]);
			fprintf(fp,"%10.6f\n", s1h[i*ND+j]);
		}
	}
	fclose(fp);
}
