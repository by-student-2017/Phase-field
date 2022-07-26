#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//#include "wingxa.h"

#define NDP 101
#define PI 3.14159
#define RR 8.3145
#define INX 600						//描画window１辺(x)のピクセルサイズ
#define INY 400						//描画window１辺の(y)ピクセルサイズ

	int nd=NDP-1;//変数sの変域(0=<s=<1)分割数
	double sh[NDP], Gch[NDP];//sとGcの配列
	double s2h[NDP], Gc2h[NDP];//sとGcの補助配列

	double Gc(double a1, double s);//Gcの関数
	void frame0(double xmin,double xmax,double dx,double ymin,double ymax,double dy);//グラフ枠描画
	void datsave();//データ保存サブルーチン
	void datload();//データ読み出しサブルーチン

//*** 自由エネルギ－描画 ***********************************************************
int main(void){

	int i, i1, i2;
	double a1, b1, c1;
	double s, s_max, Gc_max;

	double xmax, xmin, dx, ymax, ymin, dy;
	int ixmin=0, ixmax=INX; // 描画windowの幅
	int iymin=0, iymax=INY; // 描画windowの高さ
	double gx, gy, gx1, gy1, gx2, gy2;
	int igx, igy, igx1, igy1, igx2, igy2;

	xmin=0.0; xmax=1.0 ; dx=0.1;//sの描画範囲
	ymin=-1.0; ymax=1.0; dy=0.1;//Ｇｃの描画範囲

//--- デ－タの計算 -------------------------------------------------------
	printf("a1(10.0) =  "); scanf(" %lf",&a1);//aの入力

	s_max=a1/(2.0*a1+12.0);//Gcの極大値を与えるsの値
	Gc_max=a1*a1*a1*(a1+8.0)/32.0/pow((a1+6.0),3.0);//Gcの極大値
	printf("s_max, Gc_max= %f, %e \n", s_max, Gc_max);//標準入出力に数値を表示

	for(i=0;i<=nd;i++){
		s=(double)i/(double)nd; sh[i]=s;  Gch[i]=Gc(a1, s);//データ(s,Gc)の計算
	}

//--- デ－タの描画 -------------------------------------------------
	//frame0(xmin,xmax,dx,ymin,ymax,dy);//windowの描画

	//for(i=0;i<nd;i++){
	//	i1=i; i2=i+1;			//個々の隣り合う２点を、直線で連続的に結ぶ

	//	//実際の値
	//	gx1=sh[i1];  gy1=Gch[i1];  gx2=sh[i2];  gy2=Gch[i2];

	//	//スクリーン座標に変換
	//	igx1=(ixmax-ixmin)/(xmax-xmin)*(gx1-xmin)+ixmin;
	//	igy1=(iymax-iymin)/(ymax-ymin)*(gy1-ymin)+iymin;
	//	igx2=(ixmax-ixmin)/(xmax-xmin)*(gx2-xmin)+ixmin;
	//	igy2=(iymax-iymin)/(ymax-ymin)*(gy2-ymin)+iymin;
	//	glinewidth(3);//線の幅
	//	gcolor(0,0,255); gline(igx1, igy1, igx2, igy2);//隣り合う２点を、直線で結ぶ
	//}

	datsave();//データ(s,Gc)をハードディスクに保存

	datload();//データ(s,Gc)をハードディスクから読み出し

	sleep(2000);//2秒休み

	//for(i=0;i<=nd;i++){
	//	gx=s2h[i]; gy=Gc2h[i];//ハードディスクから読み出したデータ(s,Gc)
	//	//printf("s, Gc= %f,  %e \n", gx, gy);
	//	igx=(ixmax-ixmin)/(xmax-xmin)*(gx-xmin)+ixmin;
	//	igy=(iymax-iymin)/(ymax-ymin)*(gy-ymin)+iymin;
	//	glinewidth(2);//線の幅
	//	gcolor(255,0,0); gcircle(igx, igy, 4);//データ点を中心に小さな円を描画
	//}

	//if(keypress()){return 0;}
	return 1;
}

//*** 自由エネルギー関数 *****************************************************************
double Gc(double a1, double s)
{
	double b1, c1, gc;

	b1=3.0*a1+12.0;  c1=2.0*a1+12.0;
	gc=a1*s*s/2.0-b1*s*s*s/3.0+c1*pow(s,4.0)/4.0;
	return(gc);
}

//**** windowの設定 ************************************************************
//void frame0(double xmin,double xmax,double dx,double ymin,double ymax,double dy){
//	int idx, idy, i;
//	int ixmin=0, ixmax=INX; // 描画windowの幅
//	int iymin=0, iymax=INY; // 描画windowの高さ

//	gwinsize(ixmax,iymax); ginit(1); gsetorg(0,0); //描画windowの表示
//	gcolor(255,255,255); grect(0,0, ixmax,iymax);  //window内を白塗り

//	idx=ixmax*(dx/(xmax-xmin))+0.5;  idy=iymax*(dy/(ymax-ymin))+0.5;
//	gcolor(0,0,0); grectangle(ixmin, iymin, ixmax, iymax);//windowの外枠
//	for(i=0;i<ixmax;i+=idx){gcolor(0,0,0); gline(i, iymin, i, iymax);}//window内の罫線
//	for(i=0;i<iymax;i+=idy){gcolor(0,0,0); gline(ixmin, i, ixmax, i);}//window内の罫線
//}


//**** データの保存 ************************************************
void datsave()
{
	FILE		*stream;
	int 		i;

	stream = fopen("test.dat", "w");//上書きの場合
	//stream = fopen("test.dat", "a");//追記の場合
	for(i=0;i<=nd;i++){
		fprintf(stream, "%e  %e  \n", sh[i], Gch[i]);//データの保存
	}
	fprintf(stream, "\n");
	fclose(stream);
}
//**** データの読み出し ***********************************************
void datload()
{
	FILE	*datload0;
	int i;

	datload0 = fopen("test.dat", "r");
	for(i=0;i<=nd;i++){
		fscanf(datload0, "%lf  %lf  \n", &s2h[i], &Gc2h[i]);
	}
	fclose(datload0);
}

