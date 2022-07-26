#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

#include "wingxa.h"

#define DRND(x) ((double)(x)/RAND_MAX*rand()) 	//乱数関数の設定

#define ND 80								//１辺の差分分割数
#define INXY 400						//描画window１辺のピクセルサイズ

	int nd=ND, ndm=ND-1; 			//組織の差分分割数、組織の差分分割数-1
	double rr=8.3145;					//ガス定数
	double ca, temp, time1;			//合金組成、温度、計算時間
	double ch[ND][ND];			//組織内の濃度デ－タ配列

	void shokiha_c();					//初期濃度場設定サブル－チン
	void graph_c();					//濃度場表示サブル－チン
	void datsave();					//濃度デ－タ保存サブル－チン

//******* メインプログラム ******************************************************
int main(void){

	double ck[ND][ND];					//拡散ポテンシャル
	double ch2[ND][ND]; 					//組織内の濃度デ－タ予備配列
	double mu_chem, mu_str, mu_surf;		//各ポテンシャル
	int  i, j, k, l; 								//整数
	double time1max; 							//最大時間
	double amob_c;								//原子の易動度定数
	double c;										//濃度
	double cddtt;								//濃度の増分
	double kapa_c;								//濃度勾配エネルギ－定数
	double al;									//計算領域
	double b1;									//差分ブロックサイズ
	double delt;									//時間刻み
	double a0;									//Feの格子定数
	double vm0;									//モル体積
	double c11a, c12a, c44a;					//Feの弾性率
	double c11b, c12b, c44b;					//Crの弾性率
	double c11, c12, c44;						//合金の平均弾性率
	double y100;								//弾性関数
	double eta;									//格子ミスマッチ
	double Da, Db, Dab;						//自己拡散係数とその比
	double Mc, dMc;							//易動度関数とその微分
	double L0;									//原子間相互作用パラメータ

	double cip, cim, cjp, cjm; 	//差分ブロックにおいてcを中心に、その上下左右の濃度
	double ck1dev, ck2dev;  	//拡散ポテンシャルの１階微分と２階微分
	double ck0, ckip, ckim, ckjp, ckjm; //差分ブロックにてck0を中心にその上下左右の拡散ポテンシャル
	int   ip, im, jp, jm; 		//整数（i+1, i-1, j+1, j-1）
	double sumc, dca; 			// 濃度場の総和、平均組成からのずれ

//****************************************************************************

	printf("DELT(0.2)=  "); scanf(" %lf",&delt);		//時間きざみ入力
//	delt=0.2; 					//標準値

	printf("ca(0.4) = "); scanf(" %lf",&ca); 		//合金組成入力
//	ca=0.4; 					//標準値

	printf("Temp(K)(673.) = "); scanf(" %lf",&temp); 	//時効温度入力
//	temp=673.;				//標準値

	al=150.;					//計算領域の一辺の長さ(nm)
	al=al*1.0e-9;				//(m)に変換
	b1=al/nd;					//差分１ブロックのサイズ
	amob_c=1.;				//易動度（１に規格化）
	time1=0.;					//スタ－ト時間(ル－プの回数)
	time1max=100001.;		//計算打切時間(ル－プの回数)

	a0=2.8664E-10;  				//Feの格子定数
	vm0=6.02E23*a0*a0*a0/2.; 	//モル体積
	L0=(21020.8-9.31889*temp)/rr/temp; //原子間相互作用パラメ－タ（すでに無次元化済み）[式(2-2)]
	kapa_c=6.0e-15/b1/b1/rr/temp;	//濃度勾配エネルギ－係数（ で無次元化）
	eta=0.00614;						//格子ミスマッチ

	c11a=2.331e11*vm0/rr/temp; 		//Feの弾性率、 で無次元化
	c12a=1.3544e11*vm0/rr/temp;
	c44a=1.1783e11*vm0/rr/temp;
	c11b=3.5e11*vm0/rr/temp; 		//Crの弾性率
	c12b=0.678e11*vm0/rr/temp;
	c44b=1.008e11*vm0/rr/temp;

	c11=(1.-ca)*c11a+ca*c11b;			//合金の弾性率
	c12=(1.-ca)*c12a+ca*c12b;
	c44=(1.-ca)*c44a+ca*c44b;

	y100=c11+c12-2.*(c12*c12/c11);		//弾性関数Y<100>

	Da=1.0e-4*exp(-294000./rr/temp);		//Feの自己拡散係数
	Db=2.0e-5*exp(-308000./rr/temp);		//Crの自己拡散係数
	Dab=Db/Da;								//Feの自己拡散係数で規格化

//**********************************************************************************
	shokiha_c();								//初期濃度場の設定
 	gwinsize(INXY,INXY); ginit(2); gsetorg(0,0);//描画Window表示

//*** 相分解の時間発展過程の計算スタ－ト *******************************************
start: ;

//	if((((int)(time1) % 200)==0)) {datsave();} 	//200ル－プ毎に濃度場デ－タをディスクに保存
//  現在この行の前に//があり文全体をコメント文にしている。//を外せばデ－タの保存
//  が行われるが、数値デ－タは非常に大きくなる場合があるので、//を外す場合には、
//  ハ－ドディスクの空き容量やデ－タ保存量等を良く確認した後、//を外されたい。

	if((((int)(time1) % 100)==0)) {graph_c();}	//20ル－プ毎に濃度場をビットマップ画像として画面に表示
//******[拡散ポテンシャルの計算]****************************************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ip=i+1; im=i-1; jp=j+1; jm=j-1;
			if(i==ndm) {ip=0;}  //周期的境界条件
			if(i==0) {im=ndm;}
			if(j==ndm) {jp=0;}
			if(j==0) {jm=ndm;}

			c=ch[i][j];
			cip=ch[ip][j]; cim=ch[im][j];
			cjp=ch[i][jp]; cjm=ch[i][jm];

		 	mu_chem=L0*(1.-2.*c)+log(c)-log(1.-c); 		//化学拡散ポテンシャル[式(2-3)]
			mu_surf=-2.*kapa_c*(cip+cim+cjp+cjm-4.*c);	//勾配拡散ポテンシャル[式(2-5)]
			mu_str=2.*eta*eta*y100*(c-ca); 					//弾性拡散ポテンシャル[式(2-8)]

			ck[i][j]=mu_chem+mu_str+mu_surf; 				//拡散ポテンシャル

		}
	}

//******[濃度場の時間変化]**************************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ip=i+1; im=i-1; jp=j+1; jm=j-1;
			if(i==ndm) {ip=0;} 	//周期的境界条件
			if(i==0) {im=ndm;}
			if(j==ndm) {jp=0;}
			if(j==0) {jm=ndm;}

			c=ch[i][j];
			cip=ch[ip][j]; cim=ch[im][j];
			cjp=ch[i][jp]; cjm=ch[i][jm];

			ck0=ck[i][j];
			ckip=ck[ip][j]; ckim=ck[im][j];
			ckjp=ck[i][jp]; ckjm=ck[i][jm];

			ck2dev=(ckip+ckim+ckjp+ckjm)-4.*ck0; 							//拡散ポテンシャルの２階微分
			Mc=(ca+Dab*(1.-ca))*ca*(1.-ca);  										//易動度関数
			cddtt=amob_c*Mc*ck2dev; 							//非線形拡散方程式
			ch2[i][j]=ch[i][j]+cddtt*delt; 											//濃度場の時間発展

		}
	}

//*** [濃度場の収支の補正] *******************************************************
//*** 数値計算であるので、毎回 濃度場の収支の補正を行う必要がある。] *************
  sumc=0.;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			sumc+=ch2[i][j];
		}
	}
	dca=sumc/nd/nd-ca;

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ch[i][j]=ch2[i][j]-dca;
			if(ch[i][j]>=1.){ch[i][j]=0.9999;}
			if(ch[i][j]<=0.){ch[i][j]=0.0001;}
		}
	}

//******[時間増加]*************************************************
	if(keypress()){return 0;}
	time1=time1+1.;
	if (time1<time1max) {goto start;}

end:;
  return 0;
}

//************[初期濃度波]*****************************************
void shokiha_c()
{
	int i, j;
	double rnd0; 
  srand(time(NULL)); // 乱数初期化

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			rnd0=2.*DRND(1)-1.;
			ch[i][j]=ca+rnd0/50.;
		}
	}
}

//*******[相分解組織の描画]**************************************************
void graph_c()
{
	int i, j, ii, jj;
	double col;
	double c, x, xmax, xmin, y, ymax, ymin, rad0;
	int ixmin=0, iymin=0, igx, igy, irad0;
	int ixmax=INXY;
	int iymax=INXY;

	gcls(); //画面クリア
	xmin=0.; xmax=1.;
	ymin=0.; ymax=1.;

	printf("time %f\n",time1);
	rad0=1./nd/2.;
	irad0=(ixmax-ixmin)/(xmax-xmin)*rad0+1;

	for(i=0;i<=nd;i++){
		for(j=0;j<=nd;j++){
			x=1./nd*i+rad0;  igx=(ixmax-ixmin)/(xmax-xmin)*(x-xmin)+ixmin;
			y=1./nd*j+rad0;  igy=(iymax-iymin)/(ymax-ymin)*(y-ymin)+iymin;
			ii=i; jj=j;  if(i==nd){ii=0;}  if(j==nd){jj=0;}
			col=1.-ch[ii][jj];
			if(col>=1.0){col=1.0;}  if(col<=0.0){col=0.0;}
			gcolor((int)(255*col),(int)(255*col),(int)(255*col));
			grect(igx-irad0,igy-irad0,igx+irad0,igy+irad0);
		}
	}
	swapbuffers();

} 

//************[濃度波デ－タの保存]************************************
void datsave()
{
	FILE		*stream;
	int 		i, j;

	stream = fopen("test.dat", "a");		//保存ファイル名をtest.datとしている。
	fprintf(stream, "%e\n", time1);		//時間の書き込み
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			fprintf(stream, "%e  ", ch[i][j]);		//濃度場の書き込み
		}
	}
	fprintf(stream, "\n");
	fclose(stream);
}
//*********************************************************************

