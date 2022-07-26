#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

#include "wingxa.h"

#define DRND(x) ((double)(x)/RAND_MAX*rand())
#define ND 128
#define IG 7
#define INXY 400						//描画window１辺のピクセルサイズ

	int nd=ND, ndm=ND-1; 	//組織の分割数
	int nd2=ND/2;				 	//(組織の分割数)／2：フ－リエ変換内で使用
	int ig=IG;
	double PI=3.14159;		//円周率
	double rr=8.3145;			//ガス定数
	double time1;					//時間

	double s1h[ND][ND], s2h[ND][ND];		//マルテンサイト場

	double qs;					//フ－リエ変換(qs:-1)と逆フ－リエ変換(qs:1)の区別
	double xi[ND][ND], xr[ND][ND], xif[ND], xrf[ND];//フ－リエ変換の実部・虚部配列
	double s[ND],c[ND];	//sinとcosのテーブル
	int ik[ND];					//ビット反転配列

	void shokiha();				//初期波設定サブル－チン
	void graph_s1();			//結晶度場描画サブル－チン
	void table();					//sinとcosのテーブル作成サブル－チン
	void fft();						//１次元ＦＦＴ
	void rcfft();					//２次元ＦＦＴ
	void datsave();		//デ－タ保存サブル－チン
	void datinput();		//初期波読み込み用サブル－チン

//******* Main program ******************************************
int main(void)
{
	double s1, s2;			//結晶度
	double ep11h0[ND][ND], ep22h0[ND][ND];
	double ep11qrh0[ND][ND],	ep11qih0[ND][ND];		//
	double ep22qrh0[ND][ND],	ep22qih0[ND][ND];		//
	double s1k_chem, s1k_str, s1k_su[ND][ND];					//ポテンシャル
	double s2k_chem, s2k_str, s2k_su[ND][ND];					//ポテンシャル
	double c11, c12, c44; //弾性定数
	double lam0, mu0, nu0; //弾性定数
	double eta_s1[4][4], eta_s2[4][4];
	double ec11[ND][ND], ec22[ND][ND];	//拘束歪配列
	double ep11T, ep22T;
	double ep11_0, ep22_0;
	double ep11_a, ep22_a, ep12_a, ep21_a;
	double sig11_a, sig22_a;
	double Z11ep, Z12ep, Z21ep, Z22ep;
	double sum11, sum22;
	double s1ddtt, s2ddtt;						//場の時間変化量
	double el_fac;										//弾性定数規格化変数

	int   i, j, k, l, ii, jj, kk, iii, jjj;		//整数
	int   p, q, m, n;													//整数
	int   ip, im, jp, jm, Nstep;											//整数
	double al, temp, delt;											//計算領域(0～2π)、時間、時間きざみ
	double time1max;														//最大時間（計算を止める際に使用）
	double b1, vm0, atom_n;										//規格化長さ、モル体積、単位胞内の原子数
	double smob;								//結晶変態の易動度
	double nxx, nyy, nxy, alnn;										//逆空間の基本ベクトル、ノルム

	double AA0, AA1, AA2, AA3;		//化学的駆動力定数
	double a1_c, b1_c, c1_c;
	double a1_t, b1_t, c1_t;
	double kappa_s1, kappa_s2;											//結晶度勾配エネルギ－定数
	double ds_fac;

//****** reg data ****************************************

	printf("DELT(0.1)=  ");	scanf(" %lf",&delt);		//時間きざみ入力
	//delt=0.1;

	temp=500.0;
	al=250.0*1.0E-09;	//計算領域をal[m]にセット
	b1=al/nd;

	time1=-10.0;
	time1max=1.0e+07;	//計算時間の最大値

	smob=1.0;
	ds_fac=0.01;

	AA0=1000.0/rr/temp;
	AA1=10.0;
	AA2=3.0*AA1+12.0;
	AA3=2.0*AA1+12.0;

	kappa_s1=kappa_s2=5.0e-15/rr/temp/b1/b1;

	a1_c=b1_c=c1_c=3.563E-10;
	//a1_t=b1_t=3.541E-10; c1_t=0.5*7.216E-10;
	atom_n=4.0;  vm0=6.02E23*a1_c*b1_c*c1_c/atom_n;

//*** s1場の格子ミスマッチの設定 ***
	eta_s1[1][1]=0.08; eta_s1[2][2]=-0.04;
	eta_s1[3][3]=0.;
	eta_s1[1][2]=eta_s1[2][1]=eta_s1[1][3]=eta_s1[3][1]=eta_s1[2][3]=eta_s1[3][2]=0.;

//*** s2場の格子ミスマッチの設定 ***
	eta_s2[1][1]=eta_s1[2][2];
	eta_s2[2][2]=eta_s1[1][1];
	eta_s2[3][3]=0.;
	eta_s2[1][2]=eta_s2[2][1]=eta_s2[1][3]=eta_s2[3][1]=eta_s2[2][3]=eta_s2[3][2]=0.;

//***** Niの弾性定数 ****************************
	el_fac=1.0E+11*vm0/rr/temp;
  c11=2.508*el_fac;
  c44=1.235*el_fac;
  c12=1.500*el_fac;
  //c12=c11-2.0*c44;
	lam0=c12;		mu0=c44;
	nu0=lam0/2.0/(lam0+mu0);

printf("nu0= %f  \n", nu0);

//*** 外力の設定 ***
 	//sig22_a=0.;
 	sig22_a=-1000.*1.0e+06*vm0/rr/temp; //30MPa
	ep11_a=-lam0/4.0/mu0/(lam0+mu0)*sig22_a;
	ep22_a=(lam0+2.0*mu0)/4.0/mu0/(lam0+mu0)*sig22_a;
	//ep11_a=-0.5*lam0/mu0/(3.0*lam0+2.0*mu0)*sig22_a;
	//ep22_a=(lam0+mu0)/mu0/(3.0*lam0+2.0*mu0)*sig22_a;
	ep12_a=ep21_a=0.0;

//*** sinおよびcosテ－ブルおよび初期波の設定 ***************


	table();
	//shokiha();
	datinput();

 	gwinsize(INXY,INXY); ginit(2); gsetorg(0,0);//描画Window表示

//**** main start ******************************
start: ;

	//if(time1<=100.){Nstep=5;} else{Nstep=200;}
	//if((((int)(time1) % Nstep)==0)) {datsave();} //場デ－タの保存
	if((((int)(time1) % 10)==0)) {graph_s1();} //場デ－タの表示

//***** 界面ポテンシャル [Cahn-Hilliard]***********************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ip=i+1; im=i-1; jp=j+1; jm=j-1;
			if(i==ndm){ip=0;}  if(i==0){im=ndm;}
			if(j==ndm){jp=0;}  if(j==0){jm=ndm;}
			s1k_su[i][j]=-kappa_s1*(s1h[ip][j]+s1h[im][j]+s1h[i][jp]+s1h[i][jm]-4.0*s1h[i][j]);
			s2k_su[i][j]=-kappa_s2*(s2h[ip][j]+s2h[im][j]+s2h[i][jp]+s2h[i][jm]-4.0*s2h[i][j]);
		}
	}

//**** eigen歪場のフ－リエ変換 ep11 ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			xr[i][j]=ep11h0[i][j]=eta_s1[1][1]*s1h[i][j]+eta_s2[1][1]*s2h[i][j];
			xi[i][j]=0.;
		}
	}
	qs=-1.; rcfft();
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){ 
			ep11qrh0[i][j]=xr[i][j];
			ep11qih0[i][j]=xi[i][j];
		}
	}
	ep11qrh0[0][0]=ep11qih0[0][0]=0.;

//**** eigen歪場のフ－リエ変換 ep22 ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			xr[i][j]=ep22h0[i][j]=eta_s1[2][2]*s1h[i][j]+eta_s2[2][2]*s2h[i][j];
			xi[i][j]=0.;
		}
	}
	qs=-1.; rcfft();
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){ 
			ep22qrh0[i][j]=xr[i][j];
			ep22qih0[i][j]=xi[i][j];
		}
	}
	ep22qrh0[0][0]=ep22qih0[0][0]=0.;

//*** eigen歪場の平均値の算出 ***
	sum11=sum22=0.;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			sum11+=ep11h0[i][j];  sum22+=ep22h0[i][j];
	  }
	}
  ep11_0=sum11/nd/nd;  ep22_0=sum22/nd/nd;

//***** 拘束歪ec11の計算 *************************************
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
			xr[i][j]=Z11ep*ep11qrh0[i][j]+Z12ep*ep22qrh0[i][j];
			xi[i][j]=Z11ep*ep11qih0[i][j]+Z12ep*ep22qih0[i][j];
	 }
	}
	qs=1.; rcfft();
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ec11[i][j]=xr[i][j];
		}
	}

//*****  拘束歪ec22の計算 *****************************
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
			xr[i][j]=Z21ep*ep11qrh0[i][j]+Z22ep*ep22qrh0[i][j];
			xi[i][j]=Z21ep*ep11qih0[i][j]+Z22ep*ep22qih0[i][j];
	 }
	}
	qs=1.; rcfft();
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ec22[i][j]=xr[i][j];
		}
	}

//******  ポテンシャルの計算 ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){

			s1=s1h[i][j];  	s2=s2h[i][j];

//******  化学ポテンシャルの計算 ********************************
			s1k_chem=AA0*s1*(AA1-AA2*s1+AA3*(s1*s1+s2*s2));
			s2k_chem=AA0*s2*(AA1-AA2*s2+AA3*(s1*s1+s2*s2));

//******  弾性ポテンシャルの計算 ********************************

			ep11T=ep11h0[i][j]-ep11_0-ec11[i][j]-ep11_a;
			ep22T=ep22h0[i][j]-ep22_0-ec22[i][j]-ep22_a;

			s1k_str=ep11T*((lam0+2.0*mu0)*eta_s1[1][1]+lam0*eta_s1[2][2])
						 +ep22T*((lam0+2.0*mu0)*eta_s1[2][2]+lam0*eta_s1[1][1]);
			s2k_str=ep11T*((lam0+2.0*mu0)*eta_s2[1][1]+lam0*eta_s2[2][2])
						 +ep22T*((lam0+2.0*mu0)*eta_s2[2][2]+lam0*eta_s2[1][1]);

//******  結晶場の時間発展の計算 ********************************
			s1ddtt=-smob*(s1k_chem+s1k_su[i][j]+s1k_str);
			s2ddtt=-smob*(s2k_chem+s2k_su[i][j]+s2k_str);
			s1h[i][j]=s1h[i][j]+( s1ddtt+ds_fac*(2.0*DRND(1.)-1.0) )*delt;
			s2h[i][j]=s2h[i][j]+( s2ddtt+ds_fac*(2.0*DRND(1.)-1.0) )*delt;

//*** s1の変域は(0<=s1<=1) ***
			if(s1h[i][j]>=1.){s1h[i][j]=1.;}  if(s1h[i][j]<=0.){s1h[i][j]=0.;}
			if(s2h[i][j]>=1.){s2h[i][j]=1.;}  if(s2h[i][j]<=0.){s2h[i][j]=0.;}

		}
	}

	if(keypress()){return 0;}

	time1=time1+1.;//時間の増加
	if (time1<time1max) {goto start;}

end:;

  return 0;

}

//************ 初期波設定サブル－チン *************
void shokiha()
{
	int i, j;
  srand(time(NULL)); // 乱数初期化

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//s1h[i][j]=DRND(1.)/2.; s2h[i][j]=DRND(1.)/2.;
			if(abs(j-nd2)<(nd/40)){s1h[i][j]=DRND(1.); s2h[i][j]=DRND(1.);}
		}
	}
}

//******* phase fieldの描画 ***************************************
void graph_s1()
{
	int i, j, ii, jj;
	double col, col_R, col_G, col_B, col_RG;
	int ixmin=0, iymin=0, igx, igy, irad0;
	double c, x, xmax, xmin, y, ymax, ymin, rad0, dia0;
	int ixmax=INXY;
	int iymax=INXY;

	gcls(); //画面クリア
	xmin=0.; xmax=1.; ymin=0.; ymax=1.;

	printf("time %f\n",time1);
	dia0=1.0/nd;
	rad0=dia0/2.0;
	irad0=(ixmax-ixmin)/(xmax-xmin)*rad0+1;

	for(i=0;i<=nd;i++){
		for(j=0;j<=nd;j++){
			x=rad0+dia0*i;  igx=(ixmax-ixmin)/(xmax-xmin)*(x-xmin)+ixmin;
			y=rad0+dia0*j;  igy=(iymax-iymin)/(ymax-ymin)*(y-ymin)+iymin;
			ii=i; jj=j; if(i==nd){ii=0;} if(j==nd){jj=0;}

			col_R=s1h[ii][jj];
			col_G=s2h[ii][jj];
			col_RG=col_R+col_G;
			if(col_RG>1.){col_RG=1.;}
			col_B=1.-col_RG;

			if(col_R>=0.999){col_R=1.;} if(col_R<=0.001){col_R=0.;}
			if(col_G>=0.999){col_G=1.;} if(col_G<=0.001){col_G=0.;}
			if(col_B>=0.999){col_B=1.;} if(col_B<=0.001){col_B=0.;}

			gcolor((int)(255*col_R),(int)(255*col_G),(int)(255*col_B));
			grect(igx-irad0,igy-irad0,igx+irad0,igy+irad0);
		}
	}
	swapbuffers();
}

//******* Sin, Cos Table *************************************
void table()
{
	int it, it1, it2, mc, mn;
	double q;

	q=2.*PI/nd;
	for(it=0;it<=nd2-1;it++){
		c[it]=cos(q*it); s[it]=sin(q*it);
	}
	ik[0]=0;
	mn=nd2;
	mc=1;
	for(it1=1;it1<=ig;it1++){
		for(it2=0;it2<=mc-1;it2++){
			ik[it2+mc]=ik[it2]+mn;
		}
		mn=mn/2; mc=2*mc;
	}
}

//********** FFT **************************************
void fft()
{
	int ix, ka, kb, l2, lf, mf, n2, nf;
	double tj, tr;

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
				xrf[kb]=tr*c[ix]+tj*qs*s[ix];	xif[kb]=tj*c[ix]-tr*qs*s[ix];
			}
		}
		l2=l2*2;
	}

}

//************ RC FFT ***********************************
void rcfft()
{
	int i, ic, ir, j;

	for(ir=0;ir<=ndm;ir++){
		for(ic=0;ic<=ndm;ic++){
			xrf[ic]=xr[ir][ic];	xif[ic]=xi[ir][ic];
		}
	fft();
		for(ic=0;ic<=ndm;ic++){
			xr[ir][ic]=xrf[ik[ic]];	xi[ir][ic]=xif[ik[ic]];
		}
	}
	for(ic=0;ic<=ndm;ic++){
		for(ir=0;ir<=ndm;ir++){
			xrf[ir]=xr[ir][ic];	xif[ir]=xi[ir][ic];
		}
	fft();
		for(ir=0;ir<=ndm;ir++){
			xr[ir][ic]=xrf[ik[ir]];	xi[ir][ic]=xif[ik[ir]];
		}
	}
	if(qs>0.){return;}
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			xr[i][j]=xr[i][j]/nd/nd;	xi[i][j]=xi[i][j]/nd/nd;
		}
	}

}

//************ Data Save *******************************
void datsave()
{
	FILE		*stream;
	int 		i, j;

	stream = fopen("test.dat", "a");

	fprintf(stream, "%e\n", time1);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			fprintf(stream, "%e  %e ", s1h[i][j], s2h[i][j]);
		}
	}
	fprintf(stream, "\n");
	fclose(stream);
}

//*********** Initial Data in **************************
void datinput()
{
	FILE	*datin0;
	int i, j;
	double dami0;
	char	fname[100];

	printf("input file name = ? "); scanf("%s", fname);

	//datin0 = fopen("Mini_1.dat", "r");
	datin0 = fopen(fname, "r");
	fscanf(datin0, "%lf", &dami0);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			fscanf(datin0, "%lf  %lf  ", &s1h[i][j], &s2h[i][j]);
		}
	}
	fclose(datin0);
}

