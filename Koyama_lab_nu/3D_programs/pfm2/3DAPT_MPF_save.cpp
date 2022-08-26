//３次元計算
//粒界エネルギー密度および易動度に方位依存性なし

#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

//#include "wingxa.h"

#define DRND(x) ((double)(x)/RAND_MAX*rand())

#define ND 64
#define N 15
#define GNP 21
//#define GNP 1001
#define INXY 400					//描画window１辺のピクセルサイズ

	int nd=ND, ndm=ND-1, nm=N-1, nmm=N-2, GN=GNP-1;
	double PI=3.141592, time1;
	double RR=8.3145;
	int iout, Nstep;

	double ph[N][ND][ND][ND];		//ph[1][i][j]～ph[n][i][j]:位置[i][j]におけるpf
	int qh[N][ND][ND][ND];	 		//qh[1][i][j]～qh[n][i][j]:位置[i][j]においておける粒番号
	int n00[ND][ND][ND];			//位置[i][j]においてphが0ではない場合の個数
	int n00p[ND][ND][ND];			//位置[i][j]におけるphおよびその周辺のphが0ではない場合の個数

	double ph2[N][ND][ND][ND];
	int qh2[N][ND][ND][ND];

	double aij[GNP][GNP], wij[GNP][GNP], tij[GNP][GNP], eij[GNP][GNP];

	void ini000();					//初期濃度波設定サブル－チン
	//void graph_a();					//グラフ表示サブル－チン
	void datsave();					//濃度デ－タ保存サブル－チン
	void datsave_paraview();		//濃度デ－タ保存サブル－チン
	//void datin();					//初期濃度波読み込みサブル－チン

//******* Main program ******************************************
int main(void)
{
	int i, j, k, l, ii, jj, kk, ll, it, kkk;
	int ip, im, jp, jm, kp, km;
	int n1, n2, n3;
	int n0p;
	double time1max;
	double delt, L, b1, t, s;
	double dG, M1, W1, K1, E1;
	double mu_chem, mu_grad;
	double temp;
	double sum1, sum2, sum3, sxs;
	double pddtt;

	double gamma, delta, amobi;
	double aaa, vm0;

//****** reg data ****************************************

	printf("delt(2.0)=  "); scanf(" %lf",&delt);
//	delt=2.0;

	temp=1000.0; //K
	L=500.0;//nm
	b1=L/(double)ND*1.0e-9; //m

	vm0=7.0e-6;
  gamma=0.5*vm0/RR/temp/b1;
  delta=3.0;

  aaa=2.0/PI*sqrt(2.0*delta*gamma);
	K1=aaa*aaa;
	W1=4.0*gamma/delta;
  amobi=1.;
  M1=amobi*PI*PI/(8.0*delta);
  E1=500.0/RR/temp;
  //E1=0.0;

	time1=0.;  time1max=1.0e+20+1.0;

	for(i=1;i<=GN;i++){
		for(j=1;j<=GN;j++){
			wij[i][j]=W1;
			aij[i][j]=K1;
			tij[i][j]=M1;
			eij[i][j]=0.0;
			if( (i==GN)||(j==GN) ){eij[i][j]=E1;}
			if(i>j){eij[i][j]=-eij[i][j];}
			if(i==j){	wij[i][j]=0.0;  aij[i][j]=0.0;  tij[i][j]=0.0; eij[i][j]=0.0;}
		}
	}

//****************************************************

	ini000();
	//datin();
	//gwinsize(INXY,INXY); ginit(2); gsetorg(0,0);//描画Window表示

//**** start ******************************
Nstep = 10;
iout = -1;
start: ;
	printf("time: %f \n", time1);
	//if(time1>=200.){delt=5.0;
	//if((((int)(time1) % Nstep)==0)) {datsave();}
	if((((int)(time1) % Nstep)==0)) {datsave_paraview();}
	//if((((int)(time1) % 2)==0)) {graph_a();} 

//******  main  ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				ip=i+1; im=i-1;  jp=j+1;  jm=j-1;  kp=k+1;  km=k-1;
				if(i==ndm){ip=0;}  if(i==0){im=ndm;}
				if(j==ndm){jp=0;}  if(j==0){jm=ndm;}
				if(k==ndm){kp=0;}  if(k==0){km=ndm;}

				n0p=n00p[i][j][k];
				if(n0p>=2){
					for(n1=1; n1<=n0p; n1++){
						ii=qh[n1][i][j][k];  pddtt=0.0;
						for(n2=1; n2<=n0p; n2++){
							jj=qh[n2][i][j][k];  sum1=0.0;
							for(n3=1; n3<=n0p; n3++){
          	 		kk=qh[n3][i][j][k]; sum2=0.0;
								for(l=1;l<=n00[ip][j][k];l++){ if(qh[l][ip][j][k]==kk){sum2+=ph[l][ip][j][k];} }
								for(l=1;l<=n00[im][j][k];l++){ if(qh[l][im][j][k]==kk){sum2+=ph[l][im][j][k];} }
								for(l=1;l<=n00[i][jp][k];l++){ if(qh[l][i][jp][k]==kk){sum2+=ph[l][i][jp][k];} }
								for(l=1;l<=n00[i][jm][k];l++){ if(qh[l][i][jm][k]==kk){sum2+=ph[l][i][jm][k];} }
								for(l=1;l<=n00[i][j][kp];l++){ if(qh[l][i][j][kp]==kk){sum2+=ph[l][i][j][kp];} }
								for(l=1;l<=n00[i][j][km];l++){ if(qh[l][i][j][km]==kk){sum2+=ph[l][i][j][km];} }

           			sum1+=0.5*(aij[ii][kk]-aij[jj][kk])*(sum2-6.0*ph[n3][i][j][k])
              		                   +(wij[ii][kk]-wij[jj][kk])*ph[n3][i][j][k];
							}
          		pddtt+=-2.0*tij[ii][jj]/double(n00p[i][j][k])
            		      *(sum1-8.0/PI*eij[ii][jj]*sqrt(ph[n1][i][j][k]*ph[n2][i][j][k]));
						}

						ph2[n1][i][j][k]=ph[n1][i][j][k]+pddtt*delt;
						qh2[n1][i][j][k]=qh[n1][i][j][k];
						if(ph2[n1][i][j][k]>=1.0){ph2[n1][i][j][k]=1.0;}
						if(ph2[n1][i][j][k]<=0.0){ph2[n1][i][j][k]=0.0;}
						//if(ph2[n1][i][j][k]<=1.0e-4){ph2[n1][i][j][k]=0.0;}
					}
				}
				else{ ph2[1][i][j][k]=ph[1][i][j][k]; qh2[1][i][j][k]=qh[1][i][j][k]; }//if
			}//k
		}//j
	}//i

//--- ph2→ph, qh2→qh および ph2=0,qh2=0 -------------------------
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				for(kk=1;kk<=nm;kk++){
					ph[kk][i][j][k]=ph2[kk][i][j][k];  qh[kk][i][j][k]=qh2[kk][i][j][k];
					ph2[kk][i][j][k]=0.0;
				}
			}
		}
	}

//--- 並べ替え（降順） -------------------------
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){

				for(kk=1;kk<=nm-2;kk++){ //kk=1からkk=nmを降順(大から小)に並べ替え（kk=0は無視）
					for(l=nm-1;l>kk;l--){
						if(ph[l][i][j][k]>ph[l-1][i][j][k]){
							t=ph[l][i][j][k];  ph[l][i][j][k]=ph[l-1][i][j][k]; ph[l-1][i][j][k]=t;
							it=qh[l][i][j][k]; qh[l][i][j][k]=qh[l-1][i][j][k]; qh[l-1][i][j][k]=it;
						}
					}
				}

//--- 規格化 -------------------------
				sum1=0.0; ii=0; 
				for(kk=1;kk<=nm;kk++){ if(ph[kk][i][j][k]>0.0){ii++; sum1+=ph[kk][i][j][k];} }
				n00[i][j][k]=ii; n00p[i][j][k]=ii;
				for(kk=1;kk<=n00[i][j][k];kk++){ ph[kk][i][j][k]=ph[kk][i][j][k]/sum1; }
				//for(kk=1;kk<=n00[i][j][k];kk++){ if(sum1>=1.0){ph[kk][i][j][k]=ph[kk][i][j][k]/sum1;} }

			}
		}
	}


//--- ph=0の特別な場合(周囲にその他の方位がある場合)を追加 -------------------------
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){

				ip=i+1; im=i-1;  jp=j+1;  jm=j-1;  kp=k+1;  km=k-1;
				if(i==ndm){ip=0;}  if(i==0){im=ndm;}
				if(j==ndm){jp=0;}  if(j==0){jm=ndm;}
				if(k==ndm){kp=0;}  if(k==0){km=ndm;}

				for(kk=1;kk<=n00[ip][j][k];kk++){
					kkk=qh[kk][ip][j][k]; //位置[ip][j][k]におけるk番目の粒番号
					ii=1; for(l=1;l<=n00p[i][j][k];l++){ if(qh[l][i][j][k]==kkk){ii=0;} }
         		    //方位kkkが、位置[i][j][k]の方位に在ればii=0、無ければii=1
					if(ii==1){
						n00p[i][j][k]+=1; 
						ph[ n00p[i][j][k] ][i][j][k]=0.0; qh[ n00p[i][j][k] ][i][j][k]=kkk; 
					}
				}

				for(kk=1;kk<=n00[im][j][k];kk++){
					kkk=qh[kk][im][j][k];
					ii=1; for(l=1;l<=n00p[i][j][k];l++){ if(qh[l][i][j][k]==kkk){ii=0;} }
					if(ii==1){
						n00p[i][j][k]+=1; 
						ph[ n00p[i][j][k] ][i][j][k]=0.0; qh[ n00p[i][j][k] ][i][j][k]=kkk; 
					}
				}

				for(kk=1;kk<=n00[i][jp][k];kk++){
					kkk=qh[kk][i][jp][k];
					ii=1; for(l=1;l<=n00p[i][j][k];l++){ if(qh[l][i][j][k]==kkk){ii=0;} }
					if(ii==1){
						n00p[i][j][k]+=1; 
						ph[ n00p[i][j][k] ][i][j][k]=0.0; qh[ n00p[i][j][k] ][i][j][k]=kkk; 
					}
				}

				for(kk=1;kk<=n00[i][jm][k];kk++){
					kkk=qh[kk][i][jm][k];
					ii=1; for(l=1;l<=n00p[i][j][k];l++){ if(qh[l][i][j][k]==kkk){ii=0;} }
					if(ii==1){ 
						n00p[i][j][k]+=1; 
						ph[ n00p[i][j][k] ][i][j][k]=0.0; qh[ n00p[i][j][k] ][i][j][k]=kkk; 
					}
				}

				for(kk=1;kk<=n00[i][j][kp];kk++){
					kkk=qh[kk][i][j][kp];
					ii=1; for(l=1;l<=n00p[i][j][k];l++){ if(qh[l][i][j][k]==kkk){ii=0;} }
					if(ii==1){
						n00p[i][j][k]+=1; 
						ph[ n00p[i][j][k] ][i][j][k]=0.0; qh[ n00p[i][j][k] ][i][j][k]=kkk; 
					}
				}

				for(kk=1;kk<=n00[i][j][km];kk++){
					kkk=qh[kk][i][j][km];
					ii=1; for(l=1;l<=n00p[i][j][k];l++){ if(qh[l][i][j][k]==kkk){ii=0;} }
					if(ii==1){ 
						n00p[i][j][k]+=1; 
						ph[ n00p[i][j][k] ][i][j][k]=0.0; qh[ n00p[i][j][k] ][i][j][k]=kkk; 
					}
				}

			}
		}
	}

//******[時間増加]*************************************************
	//if(keypress()){return 0;}
	time1=time1+1.;
	if (time1<time1max) {goto start;}

	end:;
  return 0;
}

//************ Initial field *************
void ini000()
{
	int i, j, k, l, it;
	int ii, jj, kk, kkk, iGN;
	int ip, im, jp, jm, kp, km;
	double sum1, t;
	double r;
	int x1, y1, z1, r0;
	int is[ND][ND][ND];
 	//srand(time(NULL)); // 乱数初期化

//--- 初期化 -------------------------
//まず全部の粒番号をGNにする。
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				ph[1][i][j][k]=1.0;  //qh[1][i][j][k]=GN;
				n00[i][j][k]=1;      n00p[i][j][k]=1;
				is[i][j][k]=GN;
			}
		}
	}

	r0=5;//粒番号の異なる半径5の円を書く（円内に粒番号を記述）
	for(iGN=1;iGN<=GN-1;iGN++){
		x1=nd*DRND(1); y1=nd*DRND(1); z1=nd*DRND(1);
		while(is[x1][y1][z1]!=GN){x1=nd*DRND(1); y1=nd*DRND(1); z1=nd*DRND(1);}
		for(i=-r0;i<=(ndm+r0);i++){
			ii=i; if(i<0){ii=nd+i;}  if(i>ndm){ii=i-nd;}
			for(j=-r0;j<=(ndm+r0);j++){
				jj=j; if(j<0){jj=nd+j;}  if(j>ndm){jj=j-nd;}
				for(k=-r0;k<=(ndm+r0);k++){
					kk=k; if(k<0){kk=nd+k;}  if(k>ndm){kk=k-nd;}
					r=sqrt( (double(i-x1))*(double(i-x1))
								 +(double(j-y1))*(double(j-y1))
								 +(double(k-z1))*(double(k-z1)) );
					if(r<=r0){ is[ii][jj][kk]=iGN; } 
				}
			}
		}
	}

	r0=5.0; x1=nd/2; y1=nd/2; z1=nd/2;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				r=sqrt( (double(i-x1))*(double(i-x1))
                +(double(j-y1))*(double(j-y1))
                +(double(k-z1))*(double(k-z1)) );
				//if(r<=r0){ is[i][j][k]=1; } 
			}
		}
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				qh[1][i][j][k]=is[i][j][k];
			}
		}
	}

//--- ph=0の特別な場合(周囲にその他の方位がある場合)を追加 -------------------------
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){

				ip=i+1; im=i-1;  jp=j+1;  jm=j-1;  kp=k+1;  km=k-1;
				if(i==ndm){ip=0;}  if(i==0){im=ndm;}
				if(j==ndm){jp=0;}  if(j==0){jm=ndm;}
				if(k==ndm){kp=0;}  if(k==0){km=ndm;}

				for(kk=1;kk<=n00[ip][j][k];kk++){
					kkk=qh[kk][ip][j][k]; //位置[ip][j][k]におけるkk番目の粒番号
					ii=1; for(l=1;l<=n00p[i][j][k];l++){ if(qh[l][i][j][k]==kkk){ii=0;} }
         		    //方位kkkが、位置[i][j][k]の方位に在ればii=0、無ければii=1
					if(ii==1){
						n00p[i][j][k]+=1; 
						ph[ n00p[i][j][k] ][i][j][k]=0.0; qh[ n00p[i][j][k] ][i][j][k]=kkk; 
					}
				}

				for(kk=1;kk<=n00[im][j][k];kk++){
					kkk=qh[kk][im][j][k];
					ii=1; for(l=1;l<=n00p[i][j][k];l++){ if(qh[l][i][j][k]==kkk){ii=0;} }
					if(ii==1){
						n00p[i][j][k]+=1; 
						ph[ n00p[i][j][k] ][i][j][k]=0.0; qh[ n00p[i][j][k] ][i][j][k]=kkk; 
					}
				}

				for(kk=1;kk<=n00[i][jp][k];kk++){
					kkk=qh[kk][i][jp][k];
					ii=1; for(l=1;l<=n00p[i][j][k];l++){ if(qh[l][i][j][k]==kkk){ii=0;} }
					if(ii==1){
						n00p[i][j][k]+=1; 
						ph[ n00p[i][j][k] ][i][j][k]=0.0; qh[ n00p[i][j][k] ][i][j][k]=kkk; 
					}
				}

				for(kk=1;kk<=n00[i][jm][k];kk++){
					kkk=qh[kk][i][jm][k];
					ii=1; for(l=1;l<=n00p[i][j][k];l++){ if(qh[l][i][j][k]==kkk){ii=0;} }
					if(ii==1){ 
						n00p[i][j][k]+=1; 
						ph[ n00p[i][j][k] ][i][j][k]=0.0; qh[ n00p[i][j][k] ][i][j][k]=kkk; 
					}
				}

				for(kk=1;kk<=n00[i][j][kp];kk++){
					kkk=qh[kk][i][j][kp];
					ii=1; for(l=1;l<=n00p[i][j][k];l++){ if(qh[l][i][j][k]==kkk){ii=0;} }
					if(ii==1){
						n00p[i][j][k]+=1; 
						ph[ n00p[i][j][k] ][i][j][k]=0.0; qh[ n00p[i][j][k] ][i][j][k]=kkk; 
					}
				}

				for(kk=1;kk<=n00[i][j][km];kk++){
					kkk=qh[kk][i][j][km];
					ii=1; for(l=1;l<=n00p[i][j][k];l++){ if(qh[l][i][j][k]==kkk){ii=0;} }
					if(ii==1){ 
						n00p[i][j][k]+=1; 
						ph[ n00p[i][j][k] ][i][j][k]=0.0; qh[ n00p[i][j][k] ][i][j][k]=kkk; 
					}
				}

			}
		}
	}


}

//******* graph s***************************************
//void graph_a()
//{
//	int i, j, k, ii, jj, kk, kkk;
//	double col, col_R, col_G, col_B;
//	int ixmin=0, iymin=0, igx, igy, irad0;
//	int ixmax=INXY, iymax=INXY;
//	double c, x, xmax, xmin, y, ymax, ymin, rad0;

  //gcls();
//	xmin=0.; xmax=1.; ymin=0.; ymax=1.;
//	printf("time %f\n",time1);
//	rad0=1./nd/2.;  irad0=(ixmax-ixmin)/(xmax-xmin)*rad0+1;

//	for(k=0;k<=nd;k++){
//		//k=0;
//		for(i=0;i<=nd;i++){
//			for(j=0;j<=nd;j++){
//				x=1./nd*i+rad0;  igx=(ixmax-ixmin)/(xmax-xmin)*(x-xmin)+ixmin;
//				y=1./nd*j+rad0;  igy=(iymax-iymin)/(ymax-ymin)*(y-ymin)+iymin;
//				ii=i; jj=j; kk=k;
//				if(i==nd){ii=0;} if(j==nd){jj=0;} if(k==nd){kk=0;}
//					col=0.; 
//					for(kkk=1;kkk<=n00[ii][jj][kk];kkk++){ 
//						col+=pow(ph[kkk][ii][jj][kk], 4.0);
//					}
//					col_R=col_G=col_B=col;
//					//if(qh[1][ii][jj][kk]==0){col_R=1.; col_G=0.; col_B=0.;}
//					//if(qh[1][ii][jj][kk]==GN){col_R=0.; col_G=0.; col_B=1.;}
//					gcolor((int)(255*col_R),(int)(255*col_G),(int)(255*col_B));
//					grect(igx-irad0,igy-irad0,igx+irad0,igy+irad0);
//			}
//		}
//		swapbuffers();
//	}
//}

//************ Data Save *******************************
void datsave()
{
	FILE		*stream;
	int 		i, j, k, kk;
	double 	col;

	stream = fopen("test.dat", "a");
	fprintf(stream, "%e  \n", time1);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				fprintf(stream, "%d  \n", n00[i][j][k]);
				for(kk=1;kk<=n00[i][j][k];kk++){
					fprintf(stream, "%d  %e  ", qh[kk][i][j][k], ph[kk][i][j][k]);
				}
				fprintf(stream, "\n");
				//col=0.; for(kk=1;kk<=n00[i][j][k];kk++){ col+=ph[kk][i][j][k]*ph[kk][i][j][k]; }
				//fprintf(stream, "%e  ", col);
			}
		}
	}
	fprintf(stream, "\n");
	fclose(stream);
}

void datsave_paraview()
{
	FILE	*fp;
	char	fName[256];
	int 	i, j, k, kk, kk2;
	int 	flag;
	double  tph;
	
	iout = iout + 1;
	for(kk=1;kk<=N;kk++){
		sprintf(fName,"3DAPT_MPF_N%03d_result%06d.vtk",kk,iout);
		fp = fopen(fName, "w");
		fprintf(fp,"# vtk DataFile Version 3.0 \n");
		fprintf(fp,"output.vtk \n");
		fprintf(fp,"ASCII \n");
		fprintf(fp,"DATASET STRUCTURED_POINTS \n");
		fprintf(fp,"DIMENSIONS %5d %5d %5d \n",(ndm+1),(ndm+1),(ndm+1));
		fprintf(fp,"ORIGIN 0.0 0.0 0.0 \n");
		fprintf(fp,"ASPECT_RATIO 1 1 1 \n");
		fprintf(fp,"POINT_DATA %16d \n",((ndm+1)*(ndm+1)*(ndm+1)));
		fprintf(fp,"SCALARS phase_field float \n");
		fprintf(fp,"LOOKUP_TABLE default \n");
		for(k=0;k<=ndm;k++){
			for(j=0;j<=ndm;j++){
				for(i=0;i<=ndm;i++){
					flag=0;
					for(kk2=1;kk2<=n00[i][j][k];kk2++){
						if(qh[kk2][i][j][k]==kk){
							fprintf(fp,"%10.6f\n", ph[kk2][i][j][k]);
							flag=1;
						}
					}
					if(flag==0){fprintf(fp,"%10.6f\n", (0.0));}
				}
			}
		}
		fclose(fp);
	}
	for(kk=0;kk<=0;kk++){
		sprintf(fName,"3DAPT_MPF_N%03d_result%06d.vtk",kk,iout);
		fp = fopen(fName, "w");
		fprintf(fp,"# vtk DataFile Version 3.0 \n");
		fprintf(fp,"output.vtk \n");
		fprintf(fp,"ASCII \n");
		fprintf(fp,"DATASET STRUCTURED_POINTS \n");
		fprintf(fp,"DIMENSIONS %5d %5d %5d \n",(ndm+1),(ndm+1),(ndm+1));
		fprintf(fp,"ORIGIN 0.0 0.0 0.0 \n");
		fprintf(fp,"ASPECT_RATIO 1 1 1 \n");
		fprintf(fp,"POINT_DATA %16d \n",((ndm+1)*(ndm+1)*(ndm+1)));
		fprintf(fp,"SCALARS phase_field float \n");
		fprintf(fp,"LOOKUP_TABLE default \n");
		for(k=0;k<=ndm;k++){
			for(j=0;j<=ndm;j++){
				for(i=0;i<=ndm;i++){
					tph=0.0;
					for(kk2=1;kk2<=n00[i][j][k];kk2++){
						tph += ph[kk2][i][j][k]*qh[kk2][i][j][k];
					}
					fprintf(fp,"%10.6f\n", tph);
				}
			}
		}
		fclose(fp);
	}
}
//*********** Initial Data in **************************
//void datin()
//{
//	FILE		*datin0;
//	int 		i, j, k, kk;

//	datin0 = fopen("test.ini", "r");

//	for(i=0;i<=ndm;i++){
//		for(j=0;j<=ndm;j++){
//			for(k=0;k<=ndm;k++){
//				for(kk=1;kk<=nm;kk++){
//					fscanf(datin0, "%d  %lf  ", &qh[kk][i][j][k], &ph[kk][i][j][k]);
//				}
//			}
//		}
//	}
//	fclose(datin0);
//}

