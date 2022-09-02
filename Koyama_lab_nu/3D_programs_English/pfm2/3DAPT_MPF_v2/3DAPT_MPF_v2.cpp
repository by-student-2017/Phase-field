//3D calculation
//No orientation dependence of grain boundary energy density and mobility

#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

//#include "wingxa.h"

#define DRND(x) ((double)(x)/RAND_MAX*rand())

//#define ND 64
//#define N 15
//#define GNP 21
//#define GNP 1001
//#define INXY 400					//Pixel size of one side of drawing window

	//int nd=ND, ndm=ND-1, nm=N-1, nmm=N-2, GN=GNP-1;
	double PI=3.141592, time1;
	double RR=8.3145;
	int iout, Nstep;

	//double ph[N][ND][ND][ND];		//ph[1][i][j][K] to ph[n][i][j][K]: pf at position [i][j][K]
	//int qh[N][ND][ND][ND];	 	//qh[1][i][j][K] to qh[n][i][j][K]: grain number at position [i][j][K]
	//int n00[ND][ND][ND];			//number of cases where ph is not 0 at position [i][j][K]
	//int n00p[ND][ND][ND];			//Number of cases where ph at position [i][j][K] and surrounding ph is not 0

	//double ph2[N][ND][ND][ND];
	//int qh2[N][ND][ND][ND];

	//double aij[GNP][GNP], wij[GNP][GNP], tij[GNP][GNP], eij[GNP][GNP];

	void ini000(double *ph, int *qh, int *n00, int *n00p, int N, int ND, int GNP);	//Initial concentration wave setting subroutine
	//void graph_a();					//graph display subroutine
			 void datsave(double *ph, int *qh, int *n00, int N, int ND);	//Concentration data save subroutine
	void datsave_paraview(double *ph, int *qh, int *n00, int N, int ND);	//Concentration data save subroutine
	//void datin();					//Initial concentration wave reading subroutine

//******* Main program ******************************************
int main(void)
{
    int ND, N, GNP;
	int nd, ndm, nm, nmm, GN;
	
	int i, j, k, l, ii, jj, kk, ll, it, kkk;
	int ip, im, jp, jm, kp, km;
	int n1, n2, n3;
	int n0p;
	double time1max;
	double delt, L, b1, t, s;
	double dG, M1, M0, W1, W0, K1, K0, E1, E0;
	double mu_chem, mu_grad;
	double temp;
	double sum1, sum2, sum3, sxs;
	double pddtt;

	double gamma, gamma0, delta, amobi;
	double aaa, vm0;

//****** reg data ****************************************
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
	ND      = int(data[0]);
	N       = int(data[1]);
	GNP     = int(data[2]);
	delt    = data[3];	// timestep
	L       = data[4];	// [nm]
	time1max= int(data[5]);
	temp    = data[6];	// [K]
	vm0     = data[7];	// gamma=gamma0*vm0/RR/temp/b1;
	gamma0  = data[8];	// gamma=gamma0*vm0/RR/temp/b1;
	delta   = data[9];	// Not delt !!! For, aaa, W1, M1
	K0      = data[10];	// K1=K0*delta*gamma/PI/PI;
	W0      = data[11];	// W1=W0*gamma/delta;
	amobi   = data[12];	// M1=amobi*PI*PI/(M0*delta);
	M0      = data[13];	// M1=amobi*PI*PI/(M0*delta);
	E0      = data[14];	// E1=E0/RR/temp;
	Nstep   = int(data[15]);
	printf("---------------------------------\n");
	//
	nd=ND;
	ndm=ND-1;
	nm=N-1;
	nmm=N-2;
	GN=GNP-1;
	//
	//double ch[(nstep+1)*nstep]; // [(nstep+1)*nstep]=[nstep*nstep+nstep], a[i][j]= a[i*n+j] for a[][n]
	//a[z * ySize * xSize + y * xSize + x], a[i][j][k]=a[i*n*m + j*n + k]
	double *ph   = (double *)malloc(sizeof(double)*( N*ND*ND*ND + ND*ND*ND + ND*ND + ND ));	//ph[1][i][j][K] to ph[n][i][j][K]: pf at position [i][j][K]
	int    *qh   =       (int *)malloc(sizeof(int)*( N*ND*ND*ND + ND*ND*ND + ND*ND + ND ));	//qh[1][i][j][K] to qh[n][i][j][K]: grain number at position [i][j][K]
	int    *n00  =       (int *)malloc(sizeof(int)*( ND*ND*ND + ND*ND + ND ));	//number of cases where ph is not 0 at position [i][j][K]
	int    *n00p =       (int *)malloc(sizeof(int)*( ND*ND*ND + ND*ND + ND ));	//Number of cases where ph at position [i][j][K] and surrounding ph is not 0
	// n00[i*ND*ND+j*ND+k], n00p[i*ND*ND+j*ND+k]
	//
	double *ph2  = (double *)malloc(sizeof(double)*( N*ND*ND*ND + ND*ND*ND + ND*ND + ND ));	// ph2[kk*ND*ND*ND+i*ND*ND+j*ND+k];
	int    *qh2  =       (int *)malloc(sizeof(int)*( N*ND*ND*ND + ND*ND*ND + ND*ND + ND ));	// qh2[kk*ND*ND*ND+i*ND*ND+j*ND+k];
	//
	double *aij  = (double *)malloc(sizeof(double)*( GNP*GNP + GNP ));	// aij[i][j]= aij[i*GNP+j]
	double *wij  = (double *)malloc(sizeof(double)*( GNP*GNP + GNP ));	// wij[i][j]= wij[i*GNP+j]
	double *tij  = (double *)malloc(sizeof(double)*( GNP*GNP + GNP ));	// tij[i][j]= tij[i*GNP+j]
	double *eij  = (double *)malloc(sizeof(double)*( GNP*GNP + GNP ));	// eij[i][j]= eij[i*GNP+j]
	//
	//printf("delt(2.0)=  "); scanf(" %lf",&delt);	//	delt=2.0;
	//L=500.0;	// [nm]
	b1=L/(double)ND*1.0e-9;	// [m]
	
	time1=0.0;
	//time1max=1.0e+20+1.0;
	
	//temp=1000.0;	// [K]

	//vm0=7.0e-6;
	//gamma=0.5*vm0/RR/temp/b1;
	gamma=gamma0*vm0/RR/temp/b1;
	//delta=3.0;

	//aaa=2.0/PI*sqrt(2.0*delta*gamma);
	//K1=aaa*aaa;
	K1=K0*delta*gamma/PI/PI;
	//W1=4.0*gamma/delta;
	W1=W0*gamma/delta;
	//amobi=1.0;
	//M1=amobi*PI*PI/(8.0*delta);
	M1=amobi*PI*PI/(M0*delta);
	//E1=500.0/RR/temp;
	E1=E0/RR/temp;

	for(i=1;i<=GN;i++){
		for(j=1;j<=GN;j++){
			//wij[i][j]=W1;
			//aij[i][j]=K1;
			//tij[i][j]=M1;
			//eij[i][j]=0.0;
			wij[i*GNP+j]=W1;
			aij[i*GNP+j]=K1;
			tij[i*GNP+j]=M1;
			eij[i*GNP+j]=0.0;
			//if( (i==GN)||(j==GN) ){eij[i][j]=E1;}
			//if(i>j){eij[i][j]=-eij[i][j];}
			//if(i==j){	wij[i][j]=0.0;  aij[i][j]=0.0;  tij[i][j]=0.0; eij[i][j]=0.0;}
			if( (i==GN)||(j==GN) ){eij[i*GNP+j]=E1;}
			if(i>j){eij[i*GNP+j]=-eij[i*GNP+j];}
			if(i==j){wij[i*GNP+j]=0.0; aij[i*GNP+j]=0.0; tij[i*GNP+j]=0.0; eij[i*GNP+j]=0.0;}
		}
	}

//****************************************************

	ini000(ph, qh, n00, n00p, N, ND, GNP);
	//datin();
	//gwinsize(INXY,INXY); ginit(2); gsetorg(0,0);//Drawing window display

//**** start ******************************
//Nstep = 10;
iout = -1;
start: ;
	printf("time: %f \n", time1);
	//if(time1>=200.){delt=5.0;
	if((((int)(time1) % Nstep)==0)) {iout = iout + 1;}
	if((((int)(time1) % Nstep)==0)) {datsave(ph, qh, n00, N, ND);}
	if((((int)(time1) % Nstep)==0)) {datsave_paraview(ph, qh, n00, N, ND);}
	//if((((int)(time1) % 2)==0)) {graph_a();} 

//******  main  ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				ip=i+1; im=i-1;  jp=j+1;  jm=j-1;  kp=k+1;  km=k-1;
				if(i==ndm){ip=0;}  if(i==0){im=ndm;}
				if(j==ndm){jp=0;}  if(j==0){jm=ndm;}
				if(k==ndm){kp=0;}  if(k==0){km=ndm;}

				//n0p=n00p[i][j][k];
				n0p=n00p[i*ND*ND+j*ND+k];
				if(n0p>=2){
					for(n1=1; n1<=n0p; n1++){
						//ii=qh[n1][i][j][k];
						ii=qh[n1*ND*ND*ND+i*ND*ND+j*ND+k];
						pddtt=0.0;
						for(n2=1; n2<=n0p; n2++){
							//jj=qh[n2][i][j][k];
							jj=qh[n2*ND*ND*ND+i*ND*ND+j*ND+k];
							sum1=0.0;
							for(n3=1; n3<=n0p; n3++){
								//kk=qh[n3][i][j][k]; sum2=0.0;
								//for(l=1;l<=n00[ip][j][k];l++){ if(qh[l][ip][j][k]==kk){sum2+=ph[l][ip][j][k];} }
								//for(l=1;l<=n00[im][j][k];l++){ if(qh[l][im][j][k]==kk){sum2+=ph[l][im][j][k];} }
								//for(l=1;l<=n00[i][jp][k];l++){ if(qh[l][i][jp][k]==kk){sum2+=ph[l][i][jp][k];} }
								//for(l=1;l<=n00[i][jm][k];l++){ if(qh[l][i][jm][k]==kk){sum2+=ph[l][i][jm][k];} }
								//for(l=1;l<=n00[i][j][kp];l++){ if(qh[l][i][j][kp]==kk){sum2+=ph[l][i][j][kp];} }
								//for(l=1;l<=n00[i][j][km];l++){ if(qh[l][i][j][km]==kk){sum2+=ph[l][i][j][km];} }
								kk=qh[n3*ND*ND*ND+i*ND*ND+j*ND+k]; sum2=0.0;
								for(l=1;l<=n00[ip*ND*ND+j*ND+k];l++){ if(qh[l*ND*ND*ND+ip*ND*ND+j*ND+k]==kk){sum2+=ph[l*ND*ND*ND+ip*ND*ND+j*ND+k];} }
								for(l=1;l<=n00[im*ND*ND+j*ND+k];l++){ if(qh[l*ND*ND*ND+im*ND*ND+j*ND+k]==kk){sum2+=ph[l*ND*ND*ND+im*ND*ND+j*ND+k];} }
								for(l=1;l<=n00[i*ND*ND+jp*ND+k];l++){ if(qh[l*ND*ND*ND+i*ND*ND+jp*ND+k]==kk){sum2+=ph[l*ND*ND*ND+i*ND*ND+jp*ND+k];} }
								for(l=1;l<=n00[i*ND*ND+jm*ND+k];l++){ if(qh[l*ND*ND*ND+i*ND*ND+jm*ND+k]==kk){sum2+=ph[l*ND*ND*ND+i*ND*ND+jm*ND+k];} }
								for(l=1;l<=n00[i*ND*ND+j*ND+kp];l++){ if(qh[l*ND*ND*ND+i*ND*ND+j*ND+kp]==kk){sum2+=ph[l*ND*ND*ND+i*ND*ND+j*ND+kp];} }
								for(l=1;l<=n00[i*ND*ND+j*ND+km];l++){ if(qh[l*ND*ND*ND+i*ND*ND+j*ND+km]==kk){sum2+=ph[l*ND*ND*ND+i*ND*ND+j*ND+km];} }
								
								//sum1+=0.5*(aij[ii][kk]-aij[jj][kk])*(sum2-6.0*ph[n3][i][j][k])
								//	+(wij[ii][kk]-wij[jj][kk])*ph[n3][i][j][k];
								sum1+=0.5*(aij[ii*GNP+kk]-aij[jj*GNP+kk])*(sum2-6.0*ph[n3*ND*ND*ND+i*ND*ND+j*ND+k])
									+(wij[ii*GNP+kk]-wij[jj*GNP+kk])*ph[n3*ND*ND*ND+i*ND*ND+j*ND+k];
							}
							//pddtt+=-2.0*tij[ii][jj]/double(n00p[i][j][k])
							// *(sum1-8.0/PI*eij[ii][jj]*sqrt(ph[n1][i][j][k]*ph[n2][i][j][k]));
							pddtt+=-2.0*tij[ii*GNP+jj]/double(n00p[i*ND*ND+j*ND+k])
							 *(sum1-8.0/PI*eij[ii*GNP+jj]*sqrt(ph[n1*ND*ND*ND+i*ND*ND+j*ND+k]*ph[n2*ND*ND*ND+i*ND*ND+j*ND+k]));
						}
						//ph2[n1][i][j][k]=ph[n1][i][j][k]+pddtt*delt;
						//qh2[n1][i][j][k]=qh[n1][i][j][k];
						//if(ph2[n1][i][j][k]>=1.0){ph2[n1][i][j][k]=1.0;}
						//if(ph2[n1][i][j][k]<=0.0){ph2[n1][i][j][k]=0.0;}
						//if(ph2[n1][i][j][k]<=1.0e-4){ph2[n1][i][j][k]=0.0;}
						//
						ph2[n1*ND*ND*ND+i*ND*ND+j*ND+k]=ph[n1*ND*ND*ND+i*ND*ND+j*ND+k]+pddtt*delt;
						qh2[n1*ND*ND*ND+i*ND*ND+j*ND+k]=qh[n1*ND*ND*ND+i*ND*ND+j*ND+k];
						if(ph2[n1*ND*ND*ND+i*ND*ND+j*ND+k]>=1.0){ph2[n1*ND*ND*ND+i*ND*ND+j*ND+k]=1.0;}
						if(ph2[n1*ND*ND*ND+i*ND*ND+j*ND+k]<=0.0){ph2[n1*ND*ND*ND+i*ND*ND+j*ND+k]=0.0;}
					}
				}
				//else{ ph2[1][i][j][k]=ph[1][i][j][k]; qh2[1][i][j][k]=qh[1][i][j][k]; }//if
				else{
					ph2[1*ND*ND*ND+i*ND*ND+j*ND+k]=ph[1*ND*ND*ND+i*ND*ND+j*ND+k];
					qh2[1*ND*ND*ND+i*ND*ND+j*ND+k]=qh[1*ND*ND*ND+i*ND*ND+j*ND+k];
				}//if
			}//k
		}//j
	}//i
	//
//--- ph2¨ph, qh2¨qh ‚¨‚æ‚Ñ ph2=0,qh2=0 -------------------------
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				for(kk=1;kk<=nm;kk++){
					//ph[kk][i][j][k]=ph2[kk][i][j][k];  qh[kk][i][j][k]=qh2[kk][i][j][k];
					//ph2[kk][i][j][k]=0.0;
					ph[kk*ND*ND*ND+i*ND*ND+j*ND+k]=ph2[kk*ND*ND*ND+i*ND*ND+j*ND+k];
					qh[kk*ND*ND*ND+i*ND*ND+j*ND+k]=qh2[kk*ND*ND*ND+i*ND*ND+j*ND+k];
					ph2[kk*ND*ND*ND+i*ND*ND+j*ND+k]=0.0;
				}
			}
		}
	}
	//
//--- Sort (descending order) -------------------------
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//
				for(kk=1;kk<=nm-2;kk++){ //Sort from kk=1 to kk=nm in descending order (large to small) (ignoring kk=0)
					for(l=nm-1;l>kk;l--){
						//if(ph[l][i][j][k]>ph[l-1][i][j][k]){
						//	t=ph[l][i][j][k];  ph[l][i][j][k]=ph[l-1][i][j][k]; ph[l-1][i][j][k]=t;
						//	it=qh[l][i][j][k]; qh[l][i][j][k]=qh[l-1][i][j][k]; qh[l-1][i][j][k]=it;
						//}
						if(ph[l*ND*ND*ND+i*ND*ND+j*ND+k] > ph[(l-1)*ND*ND*ND+i*ND*ND+j*ND+k]){
							t=ph[l*ND*ND*ND+i*ND*ND+j*ND+k];
							  ph[l*ND*ND*ND+i*ND*ND+j*ND+k]=ph[(l-1)*ND*ND*ND+i*ND*ND+j*ND+k];
															ph[(l-1)*ND*ND*ND+i*ND*ND+j*ND+k]=t;
							it=qh[l*ND*ND*ND+i*ND*ND+j*ND+k];
							   qh[l*ND*ND*ND+i*ND*ND+j*ND+k]=qh[(l-1)*ND*ND*ND+i*ND*ND+j*ND+k];
															 qh[(l-1)*ND*ND*ND+i*ND*ND+j*ND+k]=it;
						}
					}
				}
				//
//--- standardization -------------------------
				sum1=0.0; ii=0;
				//for(kk=1;kk<=nm;kk++){ if(ph[kk][i][j][k]>0.0){ii++; sum1+=ph[kk][i][j][k];} }
				//n00[i][j][k]=ii; n00p[i][j][k]=ii;
				//for(kk=1;kk<=n00[i][j][k];kk++){ ph[kk][i][j][k]=ph[kk][i][j][k]/sum1; }
				//for(kk=1;kk<=n00[i][j][k];kk++){ if(sum1>=1.0){ph[kk][i][j][k]=ph[kk][i][j][k]/sum1;} }
				//
				for(kk=1;kk<=nm;kk++){
					if(ph[kk*ND*ND*ND+i*ND*ND+j*ND+k]>0.0){
						ii++;
						sum1+=ph[kk*ND*ND*ND+i*ND*ND+j*ND+k];
					}
				}
				n00[i*ND*ND+j*ND+k]=ii; n00p[i*ND*ND+j*ND+k]=ii;
				for(kk=1;kk<=n00[i*ND*ND+j*ND+k];kk++){
					ph[kk*ND*ND*ND+i*ND*ND+j*ND+k]=ph[kk*ND*ND*ND+i*ND*ND+j*ND+k]/sum1;
				}
			}
		}
	}
	//
//--- Added special case of ph=0 (if there are other bearings around) -------------------------
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){

				ip=i+1; im=i-1;  jp=j+1;  jm=j-1;  kp=k+1;  km=k-1;
				if(i==ndm){ip=0;}  if(i==0){im=ndm;}
				if(j==ndm){jp=0;}  if(j==0){jm=ndm;}
				if(k==ndm){kp=0;}  if(k==0){km=ndm;}

				//for(kk=1;kk<=n00[ip][j][k];kk++){
				//	kkk=qh[kk][ip][j][k]; //kth grain number at position [ip][j][k]
				//	ii=1; for(l=1;l<=n00p[i][j][k];l++){ if(qh[l][i][j][k]==kkk){ii=0;} }
					//ii=0 if the direction kkk is in the direction of the position [i][j][k], otherwise ii=1
				//	if(ii==1){
				//		n00p[i][j][k]+=1; 
				//		ph[ n00p[i][j][k] ][i][j][k]=0.0; qh[ n00p[i][j][k] ][i][j][k]=kkk; 
				//	}
				//}
				for(kk=1;kk<=n00[ip*ND*ND+j*ND+k];kk++){
					kkk=qh[kk*ND*ND*ND+ip*ND*ND+j*ND+k]; //kth grain number at position [ip][j][k]
					ii=1; for(l=1;l<=n00p[i*ND*ND+j*ND+k];l++){ if(qh[l*ND*ND*ND+i*ND*ND+j*ND+k]==kkk){ii=0;} }
					//ii=0 if the direction kkk is in the direction of the position [i][j][k], otherwise ii=1
					if(ii==1){
						n00p[i*ND*ND+j*ND+k]+=1; 
						ph[n00p[i*ND*ND+j*ND+k]*ND*ND*ND+i*ND*ND+j*ND+k]=0.0;
						qh[n00p[i*ND*ND+j*ND+k]*ND*ND*ND+i*ND*ND+j*ND+k]=kkk; 
					}
				}

				//for(kk=1;kk<=n00[im][j][k];kk++){
				//	kkk=qh[kk][im][j][k];
				//	ii=1; for(l=1;l<=n00p[i][j][k];l++){ if(qh[l][i][j][k]==kkk){ii=0;} }
				//	if(ii==1){
				//		n00p[i][j][k]+=1; 
				//		ph[ n00p[i][j][k] ][i][j][k]=0.0; qh[ n00p[i][j][k] ][i][j][k]=kkk; 
				//	}
				//}
				for(kk=1;kk<=n00[im*ND*ND+j*ND+k];kk++){
					kkk=qh[kk*ND*ND*ND+im*ND*ND+j*ND+k];
					ii=1; for(l=1;l<=n00p[i*ND*ND+j*ND+k];l++){ if(qh[l*ND*ND*ND+i*ND*ND+j*ND+k]==kkk){ii=0;} }
					if(ii==1){
						n00p[i*ND*ND+j*ND+k]+=1; 
						ph[n00p[i*ND*ND+j*ND+k]*ND*ND*ND+i*ND*ND+j*ND+k]=0.0;
						qh[n00p[i*ND*ND+j*ND+k]*ND*ND*ND+i*ND*ND+j*ND+k]=kkk; 
					}
				}

				//for(kk=1;kk<=n00[i][jp][k];kk++){
				//	kkk=qh[kk][i][jp][k];
				//	ii=1; for(l=1;l<=n00p[i][j][k];l++){ if(qh[l][i][j][k]==kkk){ii=0;} }
				//	if(ii==1){
				//		n00p[i][j][k]+=1; 
				//		ph[ n00p[i][j][k] ][i][j][k]=0.0; qh[ n00p[i][j][k] ][i][j][k]=kkk; 
				//	}
				//}
				for(kk=1;kk<=n00[i*ND*ND+jp*ND+k];kk++){
					kkk=qh[kk*ND*ND*ND+i*ND*ND+jp*ND+k];
					ii=1; for(l=1;l<=n00p[i*ND*ND+j*ND+k];l++){ if(qh[l*ND*ND*ND+i*ND*ND+j*ND+k]==kkk){ii=0;} }
					if(ii==1){
						n00p[i*ND*ND+j*ND+k]+=1; 
						ph[n00p[i*ND*ND+j*ND+k]*ND*ND*ND+i*ND*ND+j*ND+k]=0.0;
						qh[n00p[i*ND*ND+j*ND+k]*ND*ND*ND+i*ND*ND+j*ND+k]=kkk; 
					}
				}

				//for(kk=1;kk<=n00[i][jm][k];kk++){
				//	kkk=qh[kk][i][jm][k];
				//	ii=1; for(l=1;l<=n00p[i][j][k];l++){ if(qh[l][i][j][k]==kkk){ii=0;} }
				//	if(ii==1){ 
				//		n00p[i][j][k]+=1; 
				//		ph[ n00p[i][j][k] ][i][j][k]=0.0; qh[ n00p[i][j][k] ][i][j][k]=kkk; 
				//	}
				//}
				for(kk=1;kk<=n00[i*ND*ND+jm*ND+k];kk++){
					kkk=qh[kk*ND*ND*ND+i*ND*ND+jm*ND+k];
					ii=1; for(l=1;l<=n00p[i*ND*ND+j*ND+k];l++){ if(qh[l*ND*ND*ND+i*ND*ND+j*ND+k]==kkk){ii=0;} }
					if(ii==1){ 
						n00p[i*ND*ND+j*ND+k]+=1; 
						ph[n00p[i*ND*ND+j*ND+k]*ND*ND*ND+i*ND*ND+j*ND+k]=0.0;
						qh[n00p[i*ND*ND+j*ND+k]*ND*ND*ND+i*ND*ND+j*ND+k]=kkk; 
					}
				}

				//for(kk=1;kk<=n00[i][j][kp];kk++){
				//	kkk=qh[kk][i][j][kp];
				//	ii=1; for(l=1;l<=n00p[i][j][k];l++){ if(qh[l][i][j][k]==kkk){ii=0;} }
				//	if(ii==1){
				//		n00p[i][j][k]+=1; 
				//		ph[ n00p[i][j][k] ][i][j][k]=0.0; qh[ n00p[i][j][k] ][i][j][k]=kkk; 
				//	}
				//}
				for(kk=1;kk<=n00[i*ND*ND+j*ND+kp];kk++){
					kkk=qh[kk*ND*ND*ND+i*ND*ND+j*ND+kp];
					ii=1; for(l=1;l<=n00p[i*ND*ND+j*ND+k];l++){ if(qh[l*ND*ND*ND+i*ND*ND+j*ND+k]==kkk){ii=0;} }
					if(ii==1){
						n00p[i*ND*ND+j*ND+k]+=1; 
						ph[n00p[i*ND*ND+j*ND+k]*ND*ND*ND+i*ND*ND+j*ND+k]=0.0;
						qh[n00p[i*ND*ND+j*ND+k]*ND*ND*ND+i*ND*ND+j*ND+k]=kkk; 
					}
				}

				//for(kk=1;kk<=n00[i][j][km];kk++){
				//	kkk=qh[kk][i][j][km];
				//	ii=1; for(l=1;l<=n00p[i][j][k];l++){ if(qh[l][i][j][k]==kkk){ii=0;} }
				//	if(ii==1){ 
				//		n00p[i][j][k]+=1; 
				//		ph[ n00p[i][j][k] ][i][j][k]=0.0; qh[ n00p[i][j][k] ][i][j][k]=kkk; 
				//	}
				//}
				for(kk=1;kk<=n00[i*ND*ND+j*ND+km];kk++){
					kkk=qh[kk*ND*ND*ND+i*ND*ND+j*ND+km];
					ii=1; for(l=1;l<=n00p[i*ND*ND+j*ND+k];l++){ if(qh[l*ND*ND*ND+i*ND*ND+j*ND+k]==kkk){ii=0;} }
					if(ii==1){ 
						n00p[i*ND*ND+j*ND+k]+=1; 
						ph[n00p[i*ND*ND+j*ND+k]*ND*ND*ND+i*ND*ND+j*ND+k]=0.0;
						qh[n00p[i*ND*ND+j*ND+k]*ND*ND*ND+i*ND*ND+j*ND+k]=kkk; 
					}
				}
			//
			}
		}
	}
	//
//******[time increase]*************************************************
	//if(keypress()){return 0;}
	time1=time1+1.;
	if (time1<time1max) {goto start;}
	printf("Finished \n");

end:;
	std::exit(0);
	//return 0;
}

//************ Initial field *************
void ini000(double *ph, int *qh, int *n00, int *n00p, int N, int ND, int GNP)
{
	int i, j, k, l, it;
	int ii, jj, kk, kkk, iGN;
	int ip, im, jp, jm, kp, km;
	double sum1, t;
	double r;
	int x1, y1, z1, r0;
	//int is[NDX][NDY][NDZ];
	int *is = (int *)malloc(sizeof(int)*( ND*ND*ND + ND*ND + ND ));	//is[i][j][k]=is[i*ND*ND+j*ND+k]
	int nd=ND, ndm=ND-1, GN=GNP-1;
 	//srand(time(NULL)); // random number initialization

//--- Initialization -------------------------
//First, set all grain numbers to GN.
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//ph[1][i][j][k]=1.0;  //qh[1][i][j][k]=GN;
				//n00[i][j][k]=1;      n00p[i][j][k]=1;
				//is[i][j][k]=GN;
				ph[1*ND*ND*ND+i*ND*ND+j*ND+k]=1.0;  //qh[1][i][j][k]=GN;
				n00[i*ND*ND+j*ND+k]=1;
				n00p[i*ND*ND+j*ND+k]=1;
				is[i*ND*ND+j*ND+k]=GN;
			}
		}
	}

	r0=5;//Draw a circle with a radius of 5 with different grain numbers (describe the grain number inside the circle)
	//for(iGN=1;iGN<=GN-1;iGN++){
	for(iGN=1;iGN<=N;iGN++){
		x1=nd*DRND(1); y1=nd*DRND(1); z1=nd*DRND(1);
		//while(is[x1][y1][z1]!=GN){x1=nd*DRND(1); y1=nd*DRND(1); z1=nd*DRND(1);}
		while(is[x1*ND*ND+y1*ND+z1]!=GN){x1=nd*DRND(1); y1=nd*DRND(1); z1=nd*DRND(1);}
		for(i=-r0;i<=(ndm+r0);i++){
			ii=i; if(i<0){ii=nd+i;}  if(i>ndm){ii=i-nd;}
			for(j=-r0;j<=(ndm+r0);j++){
				jj=j; if(j<0){jj=nd+j;}  if(j>ndm){jj=j-nd;}
				for(k=-r0;k<=(ndm+r0);k++){
					kk=k; if(k<0){kk=nd+k;}  if(k>ndm){kk=k-nd;}
					r=sqrt( (double(i-x1))*(double(i-x1))
						   +(double(j-y1))*(double(j-y1))
						   +(double(k-z1))*(double(k-z1)) );
					//if(r<=r0){ is[ii][jj][kk]=iGN; }
					if(r<=r0){ is[ii*ND*ND+jj*ND+kk]=iGN; } 
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
				//qh[1][i][j][k]=is[i][j][k];
				qh[1*ND*ND*ND+i*ND*ND+j*ND+k]=is[i*ND*ND+j*ND+k];
			}
		}
	}

//--- Added special case of ph=0 (if there are other bearings around) -------------------------
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){

				ip=i+1; im=i-1;  jp=j+1;  jm=j-1;  kp=k+1;  km=k-1;
				if(i==ndm){ip=0;}  if(i==0){im=ndm;}
				if(j==ndm){jp=0;}  if(j==0){jm=ndm;}
				if(k==ndm){kp=0;}  if(k==0){km=ndm;}

				//for(kk=1;kk<=n00[ip][j][k];kk++){
				//	kkk=qh[kk][ip][j][k]; //kkth grain number at position [ip][j][k]
				//	ii=1; for(l=1;l<=n00p[i][j][k];l++){ if(qh[l][i][j][k]==kkk){ii=0;} }
				//	//ii=0 if the direction kkk is in the direction of the position [i][j][k], otherwise ii=1
				//	if(ii==1){
				//		n00p[i][j][k]+=1; 
				//		ph[ n00p[i][j][k] ][i][j][k]=0.0; qh[ n00p[i][j][k] ][i][j][k]=kkk; 
				//	}
				//}
				for(kk=1;kk<=n00[ip*ND*ND+j*ND+k];kk++){
					kkk=qh[kk*ND*ND*ND+ip*ND*ND+j*ND+k]; //kkth grain number at position [ip][j][k]
					ii=1; for(l=1;l<=n00p[i*ND*ND+j*ND+k];l++){ if(qh[l*ND*ND*ND+i*ND*ND+j*ND+k]==kkk){ii=0;} }
         		    //ii=0 if the direction kkk is in the direction of the position [i][j][k], otherwise ii=1
					if(ii==1){
						n00p[i*ND*ND+j*ND+k]+=1; 
						ph[n00p[i*ND*ND+j*ND+k]*ND*ND*ND+i*ND*ND+j*ND+k]=0.0;
						qh[n00p[i*ND*ND+j*ND+k]*ND*ND*ND+i*ND*ND+j*ND+k]=kkk; 
					}
				}

				//for(kk=1;kk<=n00[im][j][k];kk++){
				//	kkk=qh[kk][im][j][k];
				//	ii=1; for(l=1;l<=n00p[i][j][k];l++){ if(qh[l][i][j][k]==kkk){ii=0;} }
				//	if(ii==1){
				//		n00p[i][j][k]+=1; 
				//		ph[ n00p[i][j][k] ][i][j][k]=0.0; qh[ n00p[i][j][k] ][i][j][k]=kkk; 
				//	}
				//}
				for(kk=1;kk<=n00[im*ND*ND+j*ND+k];kk++){
					kkk=qh[kk*ND*ND*ND+im*ND*ND+j*ND+k];
					ii=1; for(l=1;l<=n00p[i*ND*ND+j*ND+k];l++){ if(qh[l*ND*ND*ND+i*ND*ND+j*ND+k]==kkk){ii=0;} }
					if(ii==1){
						n00p[i*ND*ND+j*ND+k]+=1; 
						ph[n00p[i*ND*ND+j*ND+k]*ND*ND*ND+i*ND*ND+j*ND+k]=0.0;
						qh[n00p[i*ND*ND+j*ND+k]*ND*ND*ND+i*ND*ND+j*ND+k]=kkk; 
					}
				}

				//for(kk=1;kk<=n00[i][jp][k];kk++){
				//	kkk=qh[kk][i][jp][k];
				//	ii=1; for(l=1;l<=n00p[i][j][k];l++){ if(qh[l][i][j][k]==kkk){ii=0;} }
				//	if(ii==1){
				//		n00p[i][j][k]+=1; 
				//		ph[ n00p[i][j][k] ][i][j][k]=0.0; qh[ n00p[i][j][k] ][i][j][k]=kkk; 
				//	}
				//}
				for(kk=1;kk<=n00[i*ND*ND+jp*ND+k];kk++){
					kkk=qh[kk*ND*ND*ND+i*ND*ND+jp*ND+k];
					ii=1; for(l=1;l<=n00p[i*ND*ND+j*ND+k];l++){ if(qh[l*ND*ND*ND+i*ND*ND+j*ND+k]==kkk){ii=0;} }
					if(ii==1){
						n00p[i*ND*ND+j*ND+k]+=1; 
						ph[n00p[i*ND*ND+j*ND+k]*ND*ND*ND+i*ND*ND+j*ND+k]=0.0;
						qh[n00p[i*ND*ND+j*ND+k]*ND*ND*ND+i*ND*ND+j*ND+k]=kkk; 
					}
				}

				//for(kk=1;kk<=n00[i][jm][k];kk++){
				//	kkk=qh[kk][i][jm][k];
				//	ii=1; for(l=1;l<=n00p[i][j][k];l++){ if(qh[l][i][j][k]==kkk){ii=0;} }
				//	if(ii==1){ 
				//		n00p[i][j][k]+=1; 
				//		ph[ n00p[i][j][k] ][i][j][k]=0.0; qh[ n00p[i][j][k] ][i][j][k]=kkk; 
				//	}
				//}
				for(kk=1;kk<=n00[i*ND*ND+jm*ND+k];kk++){
					kkk=qh[kk*ND*ND*ND+i*ND*ND+jm*ND+k];
					ii=1; for(l=1;l<=n00p[i*ND*ND+j*ND+k];l++){ if(qh[l*ND*ND*ND+i*ND*ND+j*ND+k]==kkk){ii=0;} }
					if(ii==1){ 
						n00p[i*ND*ND+j*ND+k]+=1; 
						ph[n00p[i*ND*ND+j*ND+k]*ND*ND*ND+i*ND*ND+j*ND+k]=0.0;
						qh[n00p[i*ND*ND+j*ND+k]*ND*ND*ND+i*ND*ND+j*ND+k]=kkk; 
					}
				}

				//for(kk=1;kk<=n00[i][j][kp];kk++){
				//	kkk=qh[kk][i][j][kp];
				//	ii=1; for(l=1;l<=n00p[i][j][k];l++){ if(qh[l][i][j][k]==kkk){ii=0;} }
				//	if(ii==1){
				//		n00p[i][j][k]+=1; 
				//		ph[ n00p[i][j][k] ][i][j][k]=0.0; qh[ n00p[i][j][k] ][i][j][k]=kkk; 
				//	}
				//}
				for(kk=1;kk<=n00[i*ND*ND+j*ND+kp];kk++){
					kkk=qh[kk*ND*ND*ND+i*ND*ND+j*ND+kp];
					ii=1; for(l=1;l<=n00p[i*ND*ND+j*ND+k];l++){ if(qh[l*ND*ND*ND+i*ND*ND+j*ND+k]==kkk){ii=0;} }
					if(ii==1){
						n00p[i*ND*ND+j*ND+k]+=1; 
						ph[n00p[i*ND*ND+j*ND+k]*ND*ND*ND+i*ND*ND+j*ND+k]=0.0;
						qh[n00p[i*ND*ND+j*ND+k]*ND*ND*ND+i*ND*ND+j*ND+k]=kkk; 
					}
				}

				//for(kk=1;kk<=n00[i][j][km];kk++){
				//	kkk=qh[kk][i][j][km];
				//	ii=1; for(l=1;l<=n00p[i][j][k];l++){ if(qh[l][i][j][k]==kkk){ii=0;} }
				//	if(ii==1){ 
				//		n00p[i][j][k]+=1; 
				//		ph[ n00p[i][j][k] ][i][j][k]=0.0; qh[ n00p[i][j][k] ][i][j][k]=kkk; 
				//	}
				//}
				for(kk=1;kk<=n00[i*ND*ND+j*ND+km];kk++){
					kkk=qh[kk*ND*ND*ND+i*ND*ND+j*ND+km];
					ii=1; for(l=1;l<=n00p[i*ND*ND+j*ND+k];l++){ if(qh[l*ND*ND*ND+i*ND*ND+j*ND+k]==kkk){ii=0;} }
					if(ii==1){ 
						n00p[i*ND*ND+j*ND+k]+=1; 
						ph[n00p[i*ND*ND+j*ND+k]*ND*ND*ND+i*ND*ND+j*ND+k]=0.0;
						qh[n00p[i*ND*ND+j*ND+k]*ND*ND*ND+i*ND*ND+j*ND+k]=kkk; 
					}
				}
			//
			}
		}
	}


}

//************ Data Save *******************************
void datsave(double *ph, int *qh, int *n00, int N, int ND)
{
	FILE		*stream;
	char	fName[256];
	int 		i, j, k, kk;
	double 	col;
	int nd=ND, ndm=ND-1;

	sprintf(fName,"data_%06d.dat",iout);
	//stream = fopen("test.dat", "w");
	stream = fopen(fName, "w");
	fprintf(stream, "%d %d %d \n", ndm, ndm, ndm);
	fprintf(stream, "%e  \n", time1);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				col=0.0;
				for(kk=1;kk<=n00[i*ND*ND+j*ND+k];kk++){
					if(qh[kk*ND*ND*ND+i*ND*ND+j*ND+k]<=N){
						col+=ph[kk*ND*ND*ND+i*ND*ND+j*ND+k];
					}
				}
				fprintf(stream, "%e  ", col);
			}
		}
	}
	//
	fprintf(stream, "%e  \n", time1);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//fprintf(stream, "%d  \n", n00[i][j][k]);
				//for(kk=1;kk<=n00[i][j][k];kk++){
				//	fprintf(stream, "%d  %e  ", qh[kk][i][j][k], ph[kk][i][j][k]);
				//}
				fprintf(stream, "%d  \n", n00[i*ND*ND+j*ND+k]);
				for(kk=1;kk<=n00[i*ND*ND+j*ND+k];kk++){
					fprintf(stream, "%d  %e  ", qh[kk*ND*ND*ND+i*ND*ND+j*ND+k], ph[kk*ND*ND*ND+i*ND*ND+j*ND+k]);
				}
				//
				//fprintf(stream, "\n");
				//col=0.; for(kk=1;kk<=n00[i][j][k];kk++){ col+=ph[kk][i][j][k]*ph[kk][i][j][k]; }
				//fprintf(stream, "%e  ", col);
				//
				fprintf(stream, "\n");
				col=0.; for(kk=1;kk<=n00[i*ND*ND+j*ND+k];kk++){ col+=ph[kk*ND*ND*ND+i*ND*ND+j*ND+k]*ph[kk*ND*ND*ND+i*ND*ND+j*ND+k]; }
				fprintf(stream, "%e  ", col);
			}
		}
	}
	fprintf(stream, "\n");
	fclose(stream);
}

void datsave_paraview(double *ph, int *qh, int *n00, int N, int ND)
{
	FILE	*fp;
	char	fName[256];
	int 	i, j, k, kk, kk2;
	int 	flag;
	double  tph;
	int nd=ND, ndm=ND-1;
	
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
					for(kk2=1;kk2<=n00[i*ND*ND+j*ND+k];kk2++){
						if(qh[kk2*ND*ND*ND+i*ND*ND+j*ND+k]==kk){
							fprintf(fp,"%10.6f\n", ph[kk2*ND*ND*ND+i*ND*ND+j*ND+k]);
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
					for(kk2=1;kk2<=n00[i*ND*ND+j*ND+k];kk2++){
						tph += ph[kk2*ND*ND*ND+i*ND*ND+j*ND+k]*qh[kk2*ND*ND*ND+i*ND*ND+j*ND+k];
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

