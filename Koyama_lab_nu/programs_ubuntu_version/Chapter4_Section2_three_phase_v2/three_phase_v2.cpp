#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

#define DRND(x) ((double)(x)/RAND_MAX*rand())  //乱数の設定

//#define ND 100					//差分計算における計算領域一辺の分割数

	//int nd=ND;					//計算領域の一辺の差分分割数(差分ブロック数)
	//int ndm=ND-1;				//ND-1を定義
	double PI=3.141592654;		//π
	double rr=8.3145;			//ガス定数
	double temp;				//絶対温度
	double time1;				//計算カウント数(時間に比例)
	double c2a, c3a;			//平均組成(1:A成分, 2:B成分, 3:C成分)
	//double c2h[ND][ND], c3h[ND][ND];//局所組成
	int Nstep, iout;

	void ini000(double *c2h, double *c3h, int ND);	//初期濃度プロファイルの設定サブルーチン
	void datsave(double *c2h, double *c3h, int ND);	//データ保存サブルーチン
	void datsave_paraview(double *c2h, double *c3h, int ND);//データ保存サブルーチン

//******* メインプログラム ******************************************
int main(void)
{
	int ND, nd, ndm;
	
	double c1, c2, c3;							//局所濃度
	//double c2h2[ND][ND], c3h2[ND][ND];			//局所濃度場の補助行列
	double c2k_chem, c2k_su;					//局所ポテンシャル
	//double c2k[ND][ND];						//局所ポテンシャル
	double c3k_chem, c3k_su;					//局所ポテンシャル
	//double c3k[ND][ND];						//局所ポテンシャル
	double dc2a, sumc2, dc3a, sumc3, sumc23;	//濃度場の収支計算に使用している変数
	double dakd2, dakd3;						//拡散ポテンシャルの二階微分
	double c2ddtt, c3ddtt;						//濃度場の時間変動量

	int   i, j;									//整数
	int   ip, im, jp, jm;						//(i+1),(i-1),(j+1),(j-1)
	double al, b1, rtemp, delt;					//計算領域一辺の長さ、差分プロック１辺の長さ、RT、時間きざみ
	double time1max;							//計算カウント数の最大値（計算終了カウント）
	double cmob22, cmob33, cmob23, cmob32;		//易動度

	double om_12, om_23, om_13;					//相互作用パラメータ
	double om_12e, om_23e, om_13e;				//相互作用パラメータ
	double kapa_c2, kapa_c3;					//濃度勾配エネルギー係数
	double kapa_c2c, kapa_c3c;					//濃度勾配エネルギー係数

//****** 計算条件および物質定数の設定 ****************************************
	printf("---------------------------------\n");
	printf("read parameters from parameters.txt\n");
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
	c2a     = data[1];
	c3a     = data[2];
	delt    = data[3];
	temp    = data[4];
	al      = data[5];	// [nm]
	cmob22  = data[6];
	cmob33  = data[7];
	cmob23  = data[8];
	om_12e  = data[9];
	om_13e  = data[10];
	om_23e  = data[11];
	kapa_c2c= data[12];
	kapa_c3c= data[13];
	time1max= int(data[14]);
	Nstep   = int(data[15]);
	printf("---------------------------------\n");
	//
	nd=ND;					//計算領域の一辺の差分分割数(差分ブロック数)
	ndm=ND-1;				//ND-1を定義
	//
	double *c2h  = (double *)malloc(sizeof(double)*( ND*ND ));	//局所組成
	double *c3h  = (double *)malloc(sizeof(double)*( ND*ND ));	//局所組成
	double *c2h2 = (double *)malloc(sizeof(double)*( ND*ND ));	//局所濃度場の補助行列
	double *c3h2 = (double *)malloc(sizeof(double)*( ND*ND ));	//局所濃度場の補助行列
	double *c2k  = (double *)malloc(sizeof(double)*( ND*ND ));	//局所ポテンシャル
	double *c3k  = (double *)malloc(sizeof(double)*( ND*ND ));	//局所ポテンシャル
	//
	//printf("C2B(0.333) =  ");	scanf(" %lf",&c2a);//標準入出力から平均組成(B成分)を入力	//	c2a=1./3.;
	//printf("C3C(0.333) =  ");	scanf(" %lf",&c3a);//標準入出力から平均組成(C成分)を入力	//	c3a=1./3.;
	//printf("delt(0.005)=  ");	scanf(" %lf",&delt);//時間刻み	//	delt=0.005;

	//temp=900.0;					//温度（K）
	rtemp=rr*temp;				//RT
	//al=100.0*1.0E-09;			//計算領域の一辺の長さ(m)
	al=al*1.0E-09;			//計算領域の一辺の長さ(m)
	b1=al/(double)ND;			//差分プロック１辺の長さ

	//cmob22=1.0;					//易動度
	//cmob33=1.0;					//易動度
	//cmob23=cmob32=-0.5;			//易動度
	cmob32=cmob23;

	//om_12=25000./rtemp; 		//相互作用パラメータ(J/molで、RTで無次元化)
	//om_13=25000./rtemp;
	//om_23=25000./rtemp;
	om_12=om_12e/rtemp; 		//相互作用パラメータ(J/molで、RTで無次元化)
	om_13=om_13e/rtemp;
	om_23=om_23e/rtemp;

	//kapa_c2=5.0e-15/b1/b1/rtemp;//濃度勾配エネルギー係数(Jm^2/molで、RTとb1^2で無次元化)
	//kapa_c3=5.0e-15/b1/b1/rtemp;//濃度勾配エネルギー係数(Jm^2/molで、RTとb1^2で無次元化)
	kapa_c2=kapa_c2c/b1/b1/rtemp;//濃度勾配エネルギー係数(Jm^2/molで、RTとb1^2で無次元化)
	kapa_c3=kapa_c3c/b1/b1/rtemp;//濃度勾配エネルギー係数(Jm^2/molで、RTとb1^2で無次元化)

	time1=0.0;					//計算カウント数の初期値
	//time1max=1.0e05+1.0;		//計算カウント数の最大値

//*** 初期濃度場の設定と描画Window表示 *****************************************

	ini000(c2h, c3h, ND);//初期濃度場の設定

//**** シミュレーションスタート ******************************
//Nstep = 2000;
iout = -1;
start: ;

	//if((((int)(time1) % Nstep)==0)) {datsave(c2h, c3h, ND);} //一定繰返しカウント毎に濃度場を保存
	if((((int)(time1) % Nstep)==0)) {datsave_paraview(c2h, c3h, ND);} //一定繰返しカウント毎に濃度場を保存

//***** ポテンシャル場の計算 ***********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ip=i+1; im=i-1; jp=j+1; jm=j-1;
			if(i==ndm){ip=0;}  if(i==0){im=ndm;}//周期的境界条件
			if(j==ndm){jp=0;}  if(j==0){jm=ndm;}

			//c2=c2h[i][j]; c3=c3h[i][j]; c1=1.0-c2-c3;//局所濃度場
			c2=c2h[i*ND+j]; c3=c3h[i*ND+j]; c1=1.0-c2-c3;//局所濃度場

			c2k_chem=om_12*(c1-c2)-om_13*c3+om_23*c3+(log(c2)-log(c1));//化学拡散ポテンシャル
			//c2k_su=-2.*kapa_c2*(c2h[ip][j]+c2h[im][j]+c2h[i][jp]+c2h[i][jm]-4.0*c2)
			//				  -kapa_c3*(c3h[ip][j]+c3h[im][j]+c3h[i][jp]+c3h[i][jm]-4.0*c3);//勾配ポテンシャル
			c2k_su=-2.*kapa_c2*(c2h[ip*ND+j]+c2h[im*ND+j]+c2h[i*ND+jp]+c2h[i*ND+jm]-4.0*c2)
					  -kapa_c3*(c3h[ip*ND+j]+c3h[im*ND+j]+c3h[i*ND+jp]+c3h[i*ND+jm]-4.0*c3);//勾配ポテンシャル

			c3k_chem=om_13*(c1-c3)-om_12*c2+om_23*c2+(log(c3)-log(c1));//化学拡散ポテンシャル
			//c3k_su=-2.*kapa_c3*(c3h[ip][j]+c3h[im][j]+c3h[i][jp]+c3h[i][jm]-4.0*c3)
			//				  -kapa_c2*(c2h[ip][j]+c2h[im][j]+c2h[i][jp]+c2h[i][jm]-4.0*c2);//勾配ポテンシャル
			c3k_su=-2.*kapa_c3*(c3h[ip*ND+j]+c3h[im*ND+j]+c3h[i*ND+jp]+c3h[i*ND+jm]-4.0*c3)
					  -kapa_c2*(c2h[ip*ND+j]+c2h[im*ND+j]+c2h[i*ND+jp]+c2h[i*ND+jm]-4.0*c2);//勾配ポテンシャル

			//c2k[i][j]=c2k_chem+c2k_su;//拡散ポテンシャル(式(4.1))
			//c3k[i][j]=c3k_chem+c3k_su;
			c2k[i*ND+j]=c2k_chem+c2k_su;//拡散ポテンシャル(式(4.1))
			c3k[i*ND+j]=c3k_chem+c3k_su;
		}
	}

//***** 発展方程式の計算 **********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ip=i+1; im=i-1; jp=j+1; jm=j-1;
			if(i==ndm){ip=0;} if(i==0){im=ndm;}//周期的境界条件
			if(j==ndm){jp=0;} if(j==0){jm=ndm;}

			//dakd2=c2k[ip][j]+c2k[im][j]+c2k[i][jp]+c2k[i][jm]-4.0*c2k[i][j];//拡散ポテンシャルの二階微分
			//dakd3=c3k[ip][j]+c3k[im][j]+c3k[i][jp]+c3k[i][jm]-4.0*c3k[i][j];
			dakd2=c2k[ip*ND+j]+c2k[im*ND+j]+c2k[i*ND+jp]+c2k[i*ND+jm]-4.0*c2k[i*ND+j];//拡散ポテンシャルの二階微分
			dakd3=c3k[ip*ND+j]+c3k[im*ND+j]+c3k[i*ND+jp]+c3k[i*ND+jm]-4.0*c3k[i*ND+j];

			c2ddtt=cmob22*dakd2+cmob23*dakd3;//拡散方程式(式(4.2))
			c3ddtt=cmob32*dakd2+cmob33*dakd3;

			//c2h2[i][j]=c2h[i][j]+c2ddtt*delt;//濃度場の時間発展
			//c3h2[i][j]=c3h[i][j]+c3ddtt*delt;
			c2h2[i*ND+j]=c2h[i*ND+j]+c2ddtt*delt;//濃度場の時間発展
			c3h2[i*ND+j]=c3h[i*ND+j]+c3ddtt*delt;

			//if(c2h[i][j]>=1.0){c2h[i][j]=1.0-1.0e-06;}//濃度場の変域補正
			//if(c2h[i][j]<=0.0){c2h[i][j]=1.0e-06;}
			//if(c3h[i][j]>=1.0){c3h[i][j]=1.0-1.0e-06;}
			//if(c3h[i][j]<=0.0){c3h[i][j]=1.0e-06;}
			if(c2h[i*ND+j]>=1.0){c2h[i*ND+j]=1.0-1.0e-06;}//濃度場の変域補正
			if(c2h[i*ND+j]<=0.0){c2h[i*ND+j]=1.0e-06;}
			if(c3h[i*ND+j]>=1.0){c3h[i*ND+j]=1.0-1.0e-06;}
			if(c3h[i*ND+j]<=0.0){c3h[i*ND+j]=1.0e-06;}
		}
	}

//*** 濃度場の収支補正 ***********************************************
	//sumc2=0.; for(i=0;i<=ndm;i++){  for(j=0;j<=ndm;j++){ sumc2+=c2h2[i][j]; } }
	//dc2a=sumc2/ND/ND-c2a;
	//for(i=0;i<=ndm;i++){  for(j=0;j<=ndm;j++){ c2h[i][j]=c2h2[i][j]-dc2a; } }
	sumc2=0.; for(i=0;i<=ndm;i++){  for(j=0;j<=ndm;j++){ sumc2+=c2h2[i*ND+j]; } }
	dc2a=sumc2/ND/ND-c2a;
	for(i=0;i<=ndm;i++){  for(j=0;j<=ndm;j++){ c2h[i*ND+j]=c2h2[i*ND+j]-dc2a; } }

	//sumc3=0.; for(i=0;i<=ndm;i++){  for(j=0;j<=ndm;j++){ sumc3+=c3h2[i][j]; } }
	//dc3a=sumc3/ND/ND-c3a;
	//for(i=0;i<=ndm;i++){  for(j=0;j<=ndm;j++){ c3h[i][j]=c3h2[i][j]-dc3a; } }
	sumc3=0.; for(i=0;i<=ndm;i++){  for(j=0;j<=ndm;j++){ sumc3+=c3h2[i*ND+j]; } }
	dc3a=sumc3/ND/ND-c3a;
	for(i=0;i<=ndm;i++){  for(j=0;j<=ndm;j++){ c3h[i*ND+j]=c3h2[i*ND+j]-dc3a; } }

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//sumc23=c2h[i][j]+c3h[i][j];
			//if(sumc23>=1.0){
			//	c2h[i][j]=c2h[i][j]/sumc23-1.0e-06;
			//	c3h[i][j]=c3h[i][j]/sumc23-1.0e-06;
			//}
			sumc23=c2h[i*ND+j]+c3h[i*ND+j];
			if(sumc23>=1.0){
				c2h[i*ND+j]=c2h[i*ND+j]/sumc23-1.0e-06;
				c3h[i*ND+j]=c3h[i*ND+j]/sumc23-1.0e-06;
			}
		}
	}
//*********************************************************************

	//if(keypress()){return 0;}	//キー待ち状態
	time1=time1+1.0;  if(time1<time1max){goto start;}//最大カウント数に到達したかどうかの判断

	end:;
  return 0;
}


//************ 初期濃度場の設定サブルーチン *************
void ini000(double *c2h, double *c3h, int ND)
{
	int i, j, id;
 	//srand(time(NULL));//乱数の種設定
	int ndm=ND-1;

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//c2h[i][j]=c2a+0.01*(2.0*DRND(1)-1.0);//濃度場を最大±1%の乱数にて設定
			//c3h[i][j]=c3a+0.01*(2.0*DRND(1)-1.0);
			c2h[i*ND+j]=c2a+0.01*(2.0*DRND(1)-1.0);//濃度場を最大±1%の乱数にて設定
			c3h[i*ND+j]=c3a+0.01*(2.0*DRND(1)-1.0);
		}
	}

}

//************ データ保存サブルーチン *******************************
void datsave(double *c2h, double *c3h, int ND)
{
	FILE *stream;	//ストリームのポインタ設定
	int i, j;			//整数
	int ndm=ND-1;

	stream = fopen("test.dat", "a");	//書き込む先のファイルを追記方式でオープン
	fprintf(stream, "%f\n", time1);		//計算カウント数の保存
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//fprintf(stream, "%e  %e  ", c2h[i][j], c3h[i][j]);//局所濃度場の保存
			fprintf(stream, "%e  %e  ", c2h[i*ND+j], c3h[i*ND+j]);//局所濃度場の保存
		}
	}
	fprintf(stream, "\n");	//改行の書き込み
	fclose(stream);					//ファイルをクローズ
}

void datsave_paraview(double *c2h, double *c3h, int ND)
{
	FILE	*fp;
	char	fName[256];
	int 	i, j;
	int ndm=ND-1;
	
	iout = iout + 1;
	printf("sp_result%06d.vtk \n",iout);
	sprintf(fName,"sp_result%06d.vtk",iout);
	fp = fopen(fName, "w");
	fprintf(fp,"# vtk DataFile Version 3.0 \n");
	fprintf(fp,"output.vtk \n");
	fprintf(fp,"ASCII \n");
	fprintf(fp,"DATASET STRUCTURED_POINTS \n");
	fprintf(fp,"DIMENSIONS %5d %5d %5d \n",(ndm+1),(ndm+1),(1));
	fprintf(fp,"ORIGIN 0.0 0.0 0.0 \n");
	fprintf(fp,"ASPECT_RATIO 1 1 1 \n");
	fprintf(fp,"POINT_DATA %16d \n",((ndm+1)*(ndm+1)*1));
	fprintf(fp,"SCALARS concentration_A float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			//fprintf(stream, "%e  %e  ", c2h[i][j], c3h[i][j]);//局所濃度場の保存
			//fprintf(fp,"%10.6f\n", (1.0-c2h[i][j]-c3h[i][j]));
			fprintf(fp,"%10.6f\n", (1.0-c2h[i*ND+j]-c3h[i*ND+j]));
		}
	}
	fprintf(fp,"SCALARS concentration_B float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			//fprintf(stream, "%e  %e  ", c2h[i][j], c3h[i][j]);//局所濃度場の保存
			//fprintf(fp,"%10.6f\n", c2h[i][j]);
			fprintf(fp,"%10.6f\n", c2h[i*ND+j]);
		}
	}
	fprintf(fp,"SCALARS concentration_C float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			//fprintf(stream, "%e  %e  ", c2h[i][j], c3h[i][j]);//局所濃度場の保存
			//fprintf(fp,"%10.6f\n", c3h[i][j]);
			fprintf(fp,"%10.6f\n", c3h[i*ND+j]);
		}
	}
	fclose(fp);
}