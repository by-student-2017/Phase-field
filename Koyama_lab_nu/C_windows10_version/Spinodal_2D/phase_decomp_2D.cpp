#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

#include "wingxa.h"

#define DRND(x) ((double)(x)/RAND_MAX*rand()) 	//乱数関数の設定
#define ND 64								//１辺の差分分割数
#define INXY 400						//描画window１辺のピクセルサイズ

	int nd=ND, ndm=ND-1; 			//組織の差分分割数、組織の差分分割数-1
	double rr=8.3145;					//ガス定数
	double c0, temp, time1;		//合金組成、温度、計算時間
	double ch[ND][ND];				//組織内の濃度デ－タ配列

	void ini_comp();					//初期濃度場設定サブル－チン
	void graph_c();						//濃度場表示サブル－チン
	void datsave();						//濃度デ－タ保存サブル－チン

int main(void)
{
	int  i, j, k, l; 						//整数
	int ip, im, jp, jm; 				//整数（i+1, i-1, j+1, j-1）
	double ck[ND][ND];					//拡散ポテンシャル
	double ch2[ND][ND]; 				//組織内の濃度デ－タ予備配列
	double mu_chem, mu_surf;		//各ポテンシャル
	double time1max; 						//最大時間
	double amob_c;							//原子の易動度定数
	double c;										//濃度
	double cddtt;								//濃度の増分
	double kappa_c;							//濃度勾配エネルギ－定数
	double al;									//計算領域
	double b1;									//差分ブロックサイズ
	double delt;								//時間刻み
	double a0;									//格子定数
	double Dab;									//自己拡散係数とその比
	double Mc;									//易動度関数とその微分
	double L0;									//原子間相互作用パラメータ
	double c_flu;								//濃度場の揺らぎの大きさ
	double cip, cim, cjp, cjm; 	//差分ブロックにおいてcを中心に、その上下左右の濃度
	double sumc, dc0; 					//濃度場の総和、平均組成からのずれ

//****************************************************************************

	printf("delt(0.01)=  "); scanf(" %lf",&delt);		//時間きざみ入力
	//delt=0.01; 					//標準値

	//printf("Temp(K) = "); scanf(" %lf",&temp); 		//時効温度入力
	temp=1000.;						//標準値

	printf("c0 = ");	scanf(" %lf",&c0); 						//合金組成入力
	//c0=0.3; 						//標準値
	//1000Kにおけるスピノーダル組成はおよそ0.21	
	//1000Kにおけるバイノーダル組成はおよそ0.07


	al=60.0;						//計算領域の一辺の長さ(nm)
	al=al*1.0e-9;				//(m)に変換
	b1=al/nd;						//差分１ブロックのサイズ
	time1=0.0;					//スタ－ト時間(ル－プの回数)
	time1max=1.0e+06;		//計算打切時間(ル－プの回数)

	L0=25000./rr/temp; //原子間相互作用パラメ－タ（すでに無次元化済み）[式(2-2)]
	kappa_c=5.0e-15/b1/b1/rr/temp;	//濃度勾配エネルギ－係数（ で無次元化）
	Mc=c0*(1.0-c0);  								//易動度
	c_flu=0.1;											//濃度場の揺らぎ係数

//**********************************************************************************

 	gwinsize(INXY,INXY); ginit(2); gsetorg(0,0);//描画Window表示
	ini_comp();										//初期濃度場の設定

//*** 相分解の時間発展過程の計算スタ－ト *******************************************
start: ;
	//if((((int)(time1) % 200)==0)) {datsave();}	//200ル－プ毎に濃度場を保存
	if((((int)(time1) % 200)==0)) {graph_c();}		//200ル－プ毎に濃度場を描画

//******[拡散ポテンシャルの計算]****************************************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ip=i+1; im=i-1; jp=j+1; jm=j-1;
			if(i==ndm){ip=0;}  if(i==0){im=ndm;}//周期的境界条件
			if(j==ndm){jp=0;}  if(j==0){jm=ndm;}
			c=ch[i][j]; cip=ch[ip][j]; cim=ch[im][j]; cjp=ch[i][jp]; cjm=ch[i][jm];

		 	mu_chem=L0*(1.0-2.0*c)+log(c)-log(1.-c); 			//化学ポテンシャル
			mu_surf=-2.0*kappa_c*(cip+cim+cjp+cjm-4.0*c);	//勾配ポテンシャル
			ck[i][j]=mu_chem+mu_surf; 										//拡散ポテンシャル
		}
	}

//******[濃度場の時間変化]**************************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ip=i+1; im=i-1; jp=j+1; jm=j-1;
			if(i==ndm){ip=0;}  if(i==0){im=ndm;}	//周期的境界条件
			if(j==ndm){jp=0;}  if(j==0){jm=ndm;}
			cddtt=Mc*(ck[ip][j]+ck[im][j]+ck[i][jp]+ck[i][jm]-4.0*ck[i][j]); 	//非線形拡散方程式
			ch2[i][j]=ch[i][j]+(cddtt+c_flu*(2.0*DRND(1.0)-1.0))*delt; 				//濃度場の時間発展
		}
	}

//*** [濃度場の収支の補正] *******************************************************
//*** 数値計算であるので、毎回 濃度場の収支の補正を行う必要がある。] *************
  sumc=0.;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){ sumc+=ch2[i][j]; }
	}
	dc0=sumc/nd/nd-c0;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ch[i][j]=ch2[i][j]-dc0;
			if(ch[i][j]>=1.0){ch[i][j]=0.9999;}  if(ch[i][j]<=0.0){ch[i][j]=0.0001;}
		}
	}

//******[時間増加]*************************************************

	if(keypress()){return 0;}

	time1=time1+1.0;  if(time1<time1max){goto start;}

end:;
  return 0;
}

//************[初期濃度波]*****************************************
void ini_comp()
{
	int i, j;
  srand(time(NULL)); // 乱数初期化

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ch[i][j]=c0+0.01*(2.0*DRND(1)-1.0);
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

	stream = fopen("test.dat", "a");	//保存ファイル名をtest.datとしている。
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

