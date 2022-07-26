#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

#include "wingxa.h"

#define DRND(x) ((double)(x)/RAND_MAX*rand())

#define ND 30
#define INXY 400						//描画window１辺のピクセルサイズ

	int nd=ND, ndm=ND-1; 			//組織の分割、組織の分割-1
	double rr=8.3145;
	double ca, temp, time1;
	double ch[ND][ND][ND];		//組織内の濃度デ－タ配列

	void shokiha_c();		//初期濃度波設定サブル－チン
	void graph_c();			//グラフ表示サブル－チン
	void datsave();		//濃度デ－タ保存サブル－チン
	void datin();		//初期濃度波読み込みサブル－チン

//******* Main program ******************************************************
int main(void)
{

	double ck[ND][ND][ND];//拡散ポテンシャル
	double ch2[ND][ND][ND];
	double mu_chem, mu_str, mu_surf;
	int   i, j, k, l;
	double time1max; 		//最大時間
	double amob_c;			//原子の易動度定数
	double c;				//濃度
	double cddtt;			//濃度の増分
	double kapa_c;			//濃度勾配エネルギ－定数
	double om0, om1, om2;	//原子間相互作用パラメ－タ定数
	double al;				//計算領域
	double b1;				//差分ブロックサイズ
	double delt;				//時間刻み
	double a0;				//Feの格子定数
	double vm0;				//モル体積
	double c11a, c12a, c44a;	//Feの弾性率
	double c11b, c12b, c44b;	//Crの弾性率
	double c11, c12, c44;	//合金の平均弾性率
	double y100;				//弾性関数
	double eta;				//格子ミスマッチ
	double Da, Db, Dab;		//自己拡散係数
	double Mc, dMc;			//モビリティ－関数

	double c1, c2, c2_cn=0.858;
	double Tc, Tn;
	double beta_c, beta_n;
	double c_alp;
	double d_Tc, d_Tn;
	double d_beta_c, d_beta_n;
	double d_c_alp;
	double L0;
	double fs=0.285;
	double kkf, kkp, sm;
	double d_kf, d_kp;
	double d_dg_cpm_eqm, d_dg_nm, d_gc;

	double cip, cim, cjp, cjm, ckp, ckm;
	double ck1dev, ck2dev;
	double ck0, ckip, ckim, ckjp, ckjm, ckkp, ckkm;
	int   ip, im, jp, jm, kp, km;
	double sumc, dca;

//********************************************************************************

	printf("DELT(0.05)=  "); scanf(" %lf",&delt);
//	delt=0.05;

	printf("ca (0.4)= "); scanf(" %lf",&ca);
//	ca=0.4;

	printf("Temp(K) (773.)= "); scanf(" %lf",&temp);
//	temp=773.;

	al=30.;			//(nm)
	al=al*1.0e-9;	//(m)
	b1=al/nd;
	amob_c=1.;
	time1=0.;
	time1max=100001.;

	a0=2.8664E-10;
	vm0=6.02E23*a0*a0*a0/2.;

	L0=(21020.8-9.31889*temp)/rr/temp;

	kapa_c=6.0e-15/b1/b1/rr/temp;

	eta=0.00614;

	c11a=2.331e11*vm0/rr/temp;
	c12a=1.3544e11*vm0/rr/temp;
	c44a=1.1783e11*vm0/rr/temp;

	c11b=3.5e11*vm0/rr/temp;
	c12b=0.678e11*vm0/rr/temp;
	c44b=1.008e11*vm0/rr/temp;

	c11=(1.-ca)*c11a+ca*c11b;
	c12=(1.-ca)*c12a+ca*c12b;
	c44=(1.-ca)*c44a+ca*c44b;

	y100=c11+c12-2.*(c12*c12/c11);

	Da=1.0e-4*exp(-294000./rr/temp);	//Feの自己拡散係数
	Db=2.0e-5*exp(-308000./rr/temp);	//Crの自己拡散係数
	Dab=Db/Da;							//Feの自己拡散係数で規格化

	//printf(" %f, %f %f\n", gc, dg_nm, dg_cpm_eqm);

//*************************************************************************
	shokiha_c();
 	gwinsize(INXY,INXY); ginit(2); gsetorg(0,0);//描画Window表示

//**** start **************************************************************
start: ;

//	if((((int)(time1) % 200)==0)) {datsave();}
	//if(time1==300.){datsave();}
	if((((int)(time1) % 10)==0)) {graph_c();} 

//******[拡散ポテンシャルの計算]********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){

				ip=i+1; im=i-1; jp=j+1; jm=j-1; kp=k+1; km=k-1;
				if(i==ndm) {ip=0;}  if(i==0) {im=ndm;}
				if(j==ndm) {jp=0;}  if(j==0) {jm=ndm;}
				if(k==ndm) {kp=0;}  if(k==0) {km=ndm;}

				c=ch[i][j][k];
				cip=ch[ip][j][k];  cim=ch[im][j][k];
				cjp=ch[i][jp][k];  cjm=ch[i][jm][k];
				ckp=ch[i][j][kp];  ckm=ch[i][j][km];

				mu_surf=-2.*kapa_c*(cip+cim+cjp+cjm+ckp+ckm-6.*c);

				c1=1.-c ; c2=c;

			if(c<=c2_cn){
				Tc=1043.*c1-310.*c2+(1207.3-321.2*(c2-c1))*c1*c2;
				beta_c=2.216*c1-0.4*c2+0.2525*c1*c2;
				c_alp=0.96*c1+1.0*c2;
				sm=c_alp*rr*log(beta_c+1.);
				kkf=4.*(1.-fs)*sm/(1.-exp(-4.));
				kkp=8.*fs*sm;
				d_Tc=-1353.+(1207.3-321.2*(c2-c1))*(c1-c2)-642.4*c1*c2;
				d_beta_c=-2.616+0.2525*(c1-c2);
				d_c_alp=0.04;
				d_kf=4.*(1.-fs)/(1.-exp(-4.))*d_c_alp*rr*log(beta_c+1.)
					+4.*(1.-fs)/(1.-exp(-4.))*d_beta_c*c_alp*rr/(beta_c+1.);
				d_kp=8.*fs*d_c_alp*rr*log(beta_c+1.)
					+8.*fs*d_beta_c*c_alp*rr/(beta_c+1.);
				if(temp>=Tc){d_dg_cpm_eqm=-d_kp*Tc/64.*exp(8.*(1.-temp/Tc))
										   -d_Tc*kkp/64.*exp(8.*(1.-temp/Tc))
										   -d_Tc*kkp*temp/8./Tc*exp(8.*(1.-temp/Tc));}
					else    {d_dg_cpm_eqm=1./8.*d_kp*(temp-9.*Tc/8.)
										  -9.*kkp/64.*d_Tc
								  		  -0.25*d_kf*(-temp+3./4.*Tc
													    +Tc/4.*exp(-4.*(1.-temp/Tc)))
							 	 	  	  -kkf/4.*d_Tc*(0.75+0.25*exp(-4.*(1.-temp/Tc))
														-temp/Tc*exp(-4.*(1.-temp/Tc)));}
			}
			else{
				Tn=-1873.1*c1+310.*c2;
				beta_n=-2.417*c1+0.4*c2;
				c_alp=0.96*c1+1.0*c2;
				sm=c_alp*rr*log(beta_n+1.);
				kkf=4.*(1.-fs)*sm/(1.-exp(-4.));
				kkp=8.*fs*sm;
				d_Tn=2183.1;
				d_beta_n=2.817;
				d_c_alp=0.04;
				d_kf=4.*(1.-fs)/(1.-exp(-4.))*d_c_alp*rr*log(beta_n+1.)
					+4.*(1.-fs)/(1.-exp(-4.))*d_beta_n*c_alp*rr/(beta_n+1.);
				d_kp=8.*fs*d_c_alp*rr*log(beta_n+1.)
					+8.*fs*d_beta_n*c_alp*rr/(beta_n+1.);
				if(temp>=Tn){d_dg_cpm_eqm=-d_kp*Tn/64.*exp(8.*(1.-temp/Tn))
										   -d_Tn*kkp/64.*exp(8.*(1.-temp/Tn))
										   -d_Tn*kkp*temp/8./Tn*exp(8.*(1.-temp/Tn));}
					else    {d_dg_cpm_eqm=1./8.*d_kp*(temp-9.*Tn/8.)
										  -9.*kkp/64.*d_Tn
								  		  -0.25*d_kf*(-temp+3./4.*Tn
														+Tn/4.*exp(-4.*(1.-temp/Tn)))
								  	  	  -kkf/4.*d_Tn*(0.75+0.25*exp(-4.*(1.-temp/Tn))
														-temp/Tn*exp(-4.*(1.-temp/Tn)));}
			}

				d_dg_nm=L0*(c1-c2);
		    mu_chem=log(c2)-log(c1)+d_dg_nm+d_dg_cpm_eqm/rr/temp;
				mu_str=2.*eta*eta*y100*(c-ca);
				ck[i][j][k]=mu_chem+mu_str+mu_surf;

			}
		}
	}

//******[濃度場の時間変化]**************************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
			ip=i+1; im=i-1; jp=j+1; jm=j-1; kp=k+1; km=k-1;
			if(i==ndm) {ip=0;}  if(i==0) {im=ndm;}
			if(j==ndm) {jp=0;}  if(j==0) {jm=ndm;}
			if(k==ndm) {kp=0;}  if(k==0) {km=ndm;}

			cip=ch[ip][j][k];  cim=ch[im][j][k];
			cjp=ch[i][jp][k];  cjm=ch[i][jm][k];
			ckp=ch[i][j][kp];  ckm=ch[i][j][km];

			ck0=ck[i][j][k];
			ckip=ck[ip][j][k];  ckim=ck[im][j][k];
			ckjp=ck[i][jp][k];  ckjm=ck[i][jm][k];
			ckkp=ck[i][j][kp];  ckkm=ck[i][j][km];

			ck2dev=(ckip+ckim+ckjp+ckjm+ckkp+ckkm)-6.*ck0;
			Mc=(ca+Dab*(1.-ca))*ca*(1.-ca);
			cddtt=amob_c*Mc*ck2dev;
			ch2[i][j][k]=ch[i][j][k]+cddtt*delt;
			}
		}
	}

//******[濃度場の収支の補正]**********************************
	sumc=0.;
	for(i=1;i<=ndm;i++){
	    for(j=1;j<=ndm;j++){
	    	for(k=1;k<=ndm;k++){
				sumc+=ch2[i][j][k];
	    	}
	    }
	}
    dca=sumc/nd/nd/nd-ca;

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
	    	for(k=0;k<=ndm;k++){
				ch[i][j][k]=ch2[i][j][k]-dca;
				if(ch[i][j][k]>=1.){ch[i][j][k]=0.9999;}
				if(ch[i][j][k]<=0.){ch[i][j][k]=0.0001;}
	    	}
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
	int i, j, k;
	double rnd0; 
  	srand(time(NULL)); // 乱数初期化

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
				for(k=0;k<=ndm;k++){
				rnd0=2.*DRND(1)-1.;
				ch[i][j][k]=ca+0.01*rnd0;
			}
		}
	}
}

//*******[グラフ]**************************************************
void graph_c()
{
    int i, j, k, ii, jj, kk;
    double col;
    int ixmin=0, iymin=0, igx, igy, irad0;
		int ixmax=INXY;
		int iymax=INXY;
    double c, x, xmax, xmin, y, ymax, ymin, rad0;

   	//gcls(); //画面クリア
		xmin=0.; xmax=1.;
		ymin=0.; ymax=1.;

		printf("time %f\n",time1);
		rad0=1./nd/2.;
		irad0=(ixmax-ixmin)/(xmax-xmin)*rad0+1;

	for(k=0;k<=nd;k++){
		for(i=0;i<=nd;i++){
			for(j=0;j<=nd;j++){
				x=1./nd*i+rad0;
				igx=(ixmax-ixmin)/(xmax-xmin)*(x-xmin)+ixmin;
				y=1./nd*j+rad0;
				igy=(iymax-iymin)/(ymax-ymin)*(y-ymin)+iymin;
				ii=i; jj=j, kk=k;
				if(i==nd){ii=0;}  if(j==nd){jj=0;}  if(k==nd){kk=0;}
				col=1.-ch[ii][jj][kk];
				if(col>=1.){col=1.;}  if(col<=0.){col=0.;}
				gcolor((int)(255*col),(int)(255*col),(int)(255*col));
				grect(igx-irad0,igy-irad0,igx+irad0,igy+irad0);
			}
		}
		swapbuffers();
	}
} 

//************[濃度波デ－タの保存]************************************
void datsave()
{
	FILE		*stream;
	int 		i, j, k;

	stream = fopen("test.dat", "a");
	fprintf(stream, "%e\n", time1);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				fprintf(stream, "%e  ", ch[i][j][k]);
			}
		}
	}
	fprintf(stream, "\n");
	fclose(stream);
}
//***********[ファイルからの初期濃度波入力]****************************
void datin()
{
	FILE		*datin0;
	int 		i, j;

	datin0 = fopen("test.dat", "r");
	fscanf(datin0, "%lf", &time1);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			fscanf(datin0, "%lf", &ch[i][j]);
		}
	}
	fclose(datin0);
}
