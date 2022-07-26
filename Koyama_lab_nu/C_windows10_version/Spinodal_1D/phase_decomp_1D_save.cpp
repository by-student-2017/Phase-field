#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

#include "gint.fx" 			//画面初期化関数
#include "gcls.fx" 			//画面消去関数
#include "rgb.fx"				//RGB関数
#include "getscreenx.fx"		//画面Xサイズ関数
#include "getscreeny.fx"		//画面Yサイズ関数
#include "rectfill.fx"			//長方形塗潰関数
#include "line.fx"
#include "pset.fx"
#include "circle.fx"
#include "rectangle.fx"
#include "getkey.fx"			//キ－入力待ち関数

#define DRND(x) ((double)(x)/RAND_MAX*rand()) 	//乱数関数の設定

#define ND 512								//１辺の差分分割数

	int nd=ND, ndm=ND-1; 			//組織の差分分割数、組織の差分分割数-1
	double rr=8.3145;					//ガス定数
	double c0, temp, time1;			//合金組成、温度、計算時間
	double ch[ND];			//組織内の濃度デ－タ配列

	void ini_comp();					//初期濃度場設定サブル－チン
	void graph_c();					//濃度場表示サブル－チン
	void datsave();					//濃度デ－タ保存サブル－チン

//******* メインプログラム ******************************************************
main()
{
 	gint(1000,200,"graph1",rgb(1,1,1));		//描画用Window設定

	double ck[ND];					//拡散ポテンシャル
	double ch2[ND]; 					//組織内の濃度デ－タ予備配列
	double mu_chem, mu_surf;		//各ポテンシャル
	int  i, ii; 								//整数
	double time1max; 							//最大時間
	double c;										//濃度
	double cddtt;								//濃度の増分
	double kappa_c;								//濃度勾配エネルギ－定数
	double al;									//計算領域
	double b1;									//差分ブロックサイズ
	double delt;									//時間刻み
	double Mc;							//易動度関数とその微分
	double L0;									//原子間相互作用パラメータ
	double c_flu;								//濃度場の揺らぎの大きさ

	double cip, cim; 	//差分ブロックにおいてcを中心に、その上下左右の濃度
	int   ip, im; 		//整数（i+1, i-1, j+1, j-1）
	double sumc, dc0; 			// 濃度場の総和、平均組成からのずれ

//****************************************************************************

	printf("DELT(0.01)=  "); scanf(" %lf",&delt);		//時間きざみ入力
	//delt=0.01; 					//標準値

	printf("c0(0.3) = "); scanf(" %lf",&c0); 		//合金組成入力
	//c0=0.3; 					//標準値

	//printf("Temp(K) = "); scanf(" %lf",&temp); 	//時効温度入力
	temp=1000.;				//標準値

	al=500.;					//計算領域の一辺の長さ(nm)
	al=al*1.0e-9;				//(m)に変換
	b1=al/nd;					//差分１ブロックのサイズ
	time1=0.0;					//スタ－ト時間(ル－プの回数)
	time1max=1.0+5.0e+05;		//計算打切時間(ル－プの回数)

	L0=2.5e+04/rr/temp; //原子間相互作用パラメ－タ（すでに無次元化済み）
	kappa_c=5.0e-15/b1/b1/rr/temp;	//濃度勾配エネルギ－係数（ で無次元化）
	Mc=c0*(1.0-c0);  			//易動度
	c_flu=0.1;								//濃度場の揺らぎ係数

//**********************************************************************************
	ini_comp();								//初期濃度場の設定

//*** 相分解の時間発展過程の計算スタ－ト *******************************************
start: ;

	if((((int)(time1) % 1000)==0)) {datsave();}	//濃度場を保存
	if((((int)(time1) % 500)==0)){ graph_c(); }	//濃度場を表示

//******[拡散ポテンシャルの計算]****************************************************
	for(i=0;i<=ndm;i++){
		ip=i+1; im=i-1;  if(i==ndm){ip=0;}  if(i==0){im=ndm;}//周期的境界条件
		c=ch[i]; cip=ch[ip]; cim=ch[im];
	 	mu_chem=L0*(1.-2.*c)+log(c)-log(1.-c); 		//化学ポテンシャル
		mu_surf=-2.*kappa_c*(cip+cim-2.*c);	//勾配ポテンシャル
		ck[i]=mu_chem+mu_surf; 				//拡散ポテンシャル
	}

//******[濃度場の時間変化]**************************************
	for(i=0;i<=ndm;i++){
		ip=i+1; im=i-1;  if(i==ndm){ip=0;} 	if(i==0){im=ndm;}//周期的境界条件
		cddtt=Mc*(ck[ip]+ck[im]-2.0*ck[i]); 							//非線形拡散方程式
		ch2[i]=ch[i]+(cddtt+c_flu*(2.0*DRND(1.0)-1.0))*delt;
	}

//*** [濃度場の収支の補正] *******************************************************
//*** 数値計算であるので、毎回 濃度場の収支の補正を行う必要がある。] *************
  sumc=0.; for(i=0;i<=ndm;i++){ sumc=sumc+ch2[i]; }
	dc0=sumc/nd-c0;
	for(i=0;i<=ndm;i++){
		ch[i]=ch2[i]-dc0; if(ch[i]>=1.0){ch[i]=0.99999;}  if(ch[i]<=0.0){ch[i]=0.00001;}
	}

//******[時間増加]*************************************************
	time1=time1+1.0;  if(time1<time1max){goto start;}

end:;
  return EXIT_SUCCESS;
}

//************[初期濃度波]*****************************************
void ini_comp()
{
	int i, j;
  srand(time(NULL)); // 乱数初期化
	for(i=0;i<=ndm;i++){ ch[i]=c0+0.001*(2.0*DRND(1)-1.0); }
}

//*******[相分解組織の描画]**************************************************
void graph_c()
{
	int i, ii, i1, ii1, i2, ii2;
	double col;
	double c, x, xmax, xmin, y, ymax, ymin, rad0, dx, dy;
	double gx1, gy1, gx2, gy2;
	int ixmin=0, iymin=0, igx1, igy1, igx2, igy2, irad0;
	int ixmax=getscreenx()-1; // 最大x座標  Maximum-X
	int iymax=getscreeny()-1; // 最大y座標  Maximum-Y
	int idx, idy;

	xmin=0.0; xmax=1.0; dx=0.1;
	ymin=0.0; ymax=1.0; dy=0.1;
	idx=ixmax*(dx/(xmax-xmin))+0.5;
	idy=iymax*(dy/(ymax-ymin))+0.5;

	gcls(); //画面クリア

	printf("time %f\n",time1);
	rad0=1.0/nd/2.0;
	irad0=(ixmax-ixmin)/(xmax-xmin)*rad0+1;

	for(i=0;i<=nd-1;i++){
		i1=i; i2=i+1;
		gx1=1.0/nd*i1+rad0;  gx2=1.0/nd*i2+rad0;
		ii1=i1; ii2=i2; if(ii2==nd){ii2=0;}
		gy1=ch[ii1];        gy2=ch[ii2];
		igx1=(ixmax-ixmin)/(xmax-xmin)*(gx1-xmin)+ixmin;
		igy1=iymin+iymax-((iymax-iymin)/(ymax-ymin)*(gy1-ymin)+iymin);
		igx2=(ixmax-ixmin)/(xmax-xmin)*(gx2-xmin)+ixmin;
		igy2=iymin+iymax-((iymax-iymin)/(ymax-ymin)*(gy2-ymin)+iymin);
		line(igx1, igy1, igx2, igy2, 2);			//circle(igx, igy, 1, 2);
	}

	rectangle(ixmin, iymin, ixmax, iymax, 0);
	for(i=0;i<ixmax;i+=idx){line(i, iymin, i, iymax, 0);}
	for(i=0;i<iymax;i+=idy){line(ixmin, i, ixmax, i, 0);}

}

//************[濃度波デ－タの保存]************************************
void datsave()
{
	FILE		*stream;
	int 		i;

	stream = fopen("test.dat", "a");		//保存ファイル名をtest.datとしている。
	fprintf(stream, "%e \n", time1);		//時間の書き込み
	for(i=0;i<=ndm;i++){ fprintf(stream, "%e  ", ch[i]); }		//濃度場の書き込み
	fprintf(stream, "\n");
	fclose(stream);
}

