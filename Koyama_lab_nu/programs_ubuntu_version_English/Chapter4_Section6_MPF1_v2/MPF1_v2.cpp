//Program for polycrystalline grain structure formation
//Grain number is 1-nm
//grain number nm liquid phase

#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

#define DRND(x) ((double)(x)/RAND_MAX*rand())//Random number setting

//#define ND 100				//Number of divisions per side of calculation domain in difference calculation
//#define N 21					//Number of crystal orientations to consider + 1 (compared to MPF0.cpp, this value is increased)

	//int nd=ND, ndm=ND-1;		//Define the number of difference divisions
								// (number of difference blocks) on one side of the computational domain, ND-1
	//int nm=N-1, nmm=N-2;		//Define the number of crystal orientations to consider,
								// N-2 (number of crystal orientations to consider - 1)
	double PI=3.141592;			//pi
	double RR=8.3145;			//gas constant
	double time1;				//calculation count number
	//double ph[N][ND][ND], ph2[N][ND][ND];	//phase field, phase field auxiliary array
	//double aij[N][N];			//gradient energy factor
	//double wij[N][N];			//the coefficient of the penalty term
	//double tij[N][N];			//Grain boundary mobility
	//double eij[N][N];			//Driving force of grain boundary migration
	//int m00h[N][ND][ND];		//number of directions where p is not 0 at location (i,j) and its surroundings (iÅ}1,jÅ}1)
	//int n00h[ND][ND];			//Number of orientations where p is not 0 at location (i,j) and its surroundings (iÅ}1,jÅ}1)
	int Nstep, iout;

	void ini000(double *ph, int ND, int N);				//Initial field setup subroutine
	void datsave(double *ph, int ND, int N);			//data save subroutine
	void datsave_paraview(double *ph, int ND, int N);	//data save subroutine
	void datin(double *ph, int ND, int N);				//data entry subroutine

//******* main program ******************************************
int main(void)
{
	int ND, N;
	int nd, ndm, nm, nmm;
	
	int i, j, k, l, ii, jj, kk, ll, it;	//integer
	int ip, im, jp, jm;	//integer
	int n1, n2, n3;		//integer
	int n00;		//Number of orientations where p is not 0 at location (i,j) and its surroundings (iÅ}1,jÅ}1)
	//int n000;		//Number of orientations where p is not 0 at position (i,j) (n00>=n000)
	double time1max;		//Maximum calculation count (calculation end count)
	double delt, L, b1;		//Time step, length of one side of calculation area, length of one side of difference block
	double M1, M0;			//Grain boundary mobility
	double W1, W0;			//the coefficient of the penalty term
	double K1, K0;			//gradient energy factor
	double E1, E0;			//Driving force of grain boundary migration
	double temp;			//temperature
	double sum1, sum2, sum3;//Work variables for various sums
	double pddtt;			//Time change rate of phase field

	double gamma, gamma0;	//grain boundary energy density
	double delta;			//Grain boundary width (expressed by the number of differential blocks)
	double amobi;			//Grain boundary mobility
	double vm0;				//molar volume

//****** Setting calculation conditions and material constants ****************************************
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
	delt    = data[2];
	temp    = data[3];
	L       = data[4];
	vm0     = data[5];
	gamma0  = data[6];
	delta   = data[7];
	K0      = data[8];
	W0      = data[9];
	amobi   = data[10];
	M0      = data[11];
	E0      = data[12];
	time1max= int(data[13]);
	Nstep   = int(data[14]);
	printf("---------------------------------\n");
	//
	nd=ND;
	ndm=ND-1;
	nm=N-1;
	nmm=N-2;
	//
	double *ph   = (double *)malloc(sizeof(double)*( N*ND*ND + ND*ND + ND ));	//phase field
	double *ph2  = (double *)malloc(sizeof(double)*( N*ND*ND + ND*ND + ND ));	//Phase field auxiliary array
	double *aij  = (double *)malloc(sizeof(double)*( ND*ND + ND ));	//gradient energy factor
	double *wij  = (double *)malloc(sizeof(double)*( ND*ND + ND ));	//the coefficient of the penalty term
	double *tij  = (double *)malloc(sizeof(double)*( ND*ND + ND ));	//Grain boundary mobility
	double *eij  = (double *)malloc(sizeof(double)*( ND*ND + ND ));	//Driving force of grain boundary migration
	double *m00h = (double *)malloc(sizeof(double)*( N*ND*ND + ND*ND + ND ));	//number of orientations where p is not 0 at location (i,j) and its surroundings (iÅ}1,jÅ}1)
	double *n00h = (double *)malloc(sizeof(double)*( ND*ND + ND ));	//The number of orientations where p is not 0 at location (i,j) and its surroundings (iÅ}1,jÅ}1)
	//
	//printf("delt(5.0)=  "); scanf(" %lf",&delt);//Enter time step	//	delt=5.0;

	//temp=1000.0; 				//temperature (K)
	//L=2000.0;					//Side length of calculation area (nm)
	b1=L/(double)ND*1.0e-9; 	//Length of one side of difference block (m)

	//vm0=7.0e-6;				//molar volume
	//gamma=0.5*vm0/RR/temp/b1;	//Dimensionless grain boundary energy density (0.5J/m^2)
	gamma=gamma0*vm0/RR/temp/b1;//Dimensionless grain boundary energy density (0.5J/m^2)
	//delta=7.0;				//Grain boundary width (expressed by the number of differential blocks)

	//K1=8.0*delta*gamma/PI/PI;	//gradient energy factor [equation (4.40)]
	K1=K0*delta*gamma/PI/PI;	//gradient energy factor [equation (4.40)]
	//W1=4.0*gamma/delta;		//the coefficient of the penalty term [equation (4.40)]
	W1=W0*gamma/delta;			//the coefficient of the penalty term [equation (4.40)]
	//amobi=1.0;
	//M1=amobi*PI*PI/(8.0*delta);//Grain boundary mobility [equation (4.40)]
	M1=amobi*PI*PI/(M0*delta);	//Grain boundary mobility [equation (4.40)]
	//E1=50.0/RR/temp;			//Driving force of grain boundary migration
	E1=E0/RR/temp;				//Driving force of grain boundary migration

	time1=0.0;					//Initial value of calculation count
	//time1max=10000001.;		//Maximum calculation count

//*** Formula (4.32) - Setting of array (K,W,M,E) of formula (4.35) (setting for MPF0.cpp)*********
//	for(i=1;i<=nm;i++){
//		for(j=1;j<=nm;j++){
			//wij[i][j]=W1;  aij[i][j]=K1;  tij[i][j]=0.0;  eij[i][j]=0.0;
			//if( (i==nm)||(j==nm) ){eij[i][j]=E1; tij[i][j]=M1;}
			//if(i>j){eij[i][j]=-eij[i][j];}
			//if(i==j){	wij[i][j]=0.0;  aij[i][j]=0.0;  tij[i][j]=0.0; eij[i][j]=0.0;}
//			wij[i*ND+j]=W1;  aij[i*ND+j]=K1;  tij[i*ND+j]=0.0;  eij[i*ND+j]=0.0;
//			if( (i==nm)||(j==nm) ){eij[i*ND+j]=E1; tij[i*ND+j]=M1;}
//			if(i>j){eij[i*ND+j]=-eij[i*ND+j];}
//			if(i==j){ wij[i*ND+j]=0.0; aij[i*ND+j]=0.0; tij[i*ND+j]=0.0; eij[i*ND+j]=0.0;}
//		}
//	}

//*** Equation (4.36) - set array (K,W,M,E) in equation (4.39) **************************
	for(i=1;i<=nm;i++){
		for(j=1;j<=nm;j++){
			//wij[i][j]=W1;  aij[i][j]=K1;  tij[i][j]=M1;  eij[i][j]=0.0;
			//if( (i==nm)||(j==nm) ){eij[i][j]=E1;}
			//if(i>j){eij[i][j]=-eij[i][j];}
			//if(i==j){	wij[i][j]=0.0;  aij[i][j]=0.0;  tij[i][j]=0.0; eij[i][j]=0.0;}
			wij[i*ND+j]=W1;  aij[i*ND+j]=K1;  tij[i*ND+j]=M1;  eij[i*ND+j]=0.0;
			if( (i==nm)||(j==nm) ){eij[i*ND+j]=E1;}
			if(i>j){eij[i*ND+j]=-eij[i*ND+j];}
			if(i==j){ wij[i*ND+j]=0.0; aij[i*ND+j]=0.0; tij[i*ND+j]=0.0; eij[i*ND+j]=0.0;}
		}
	}

//*** Initial Field Settings and Drawing Window Display *****************************************
	for(k=1;k<=nm;k++){
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				//ph[k][i][j]=0.0;  ph2[k][i][j]=0.0;//Phase field and auxiliary array initialization
				 ph[k*ND*ND+i*ND+j]=0.0;	//Phase field initialization
				ph2[k*ND*ND+i*ND+j]=0.0;	//Auxiliary array initialization
			}
		}
	}
	ini000(ph, ND, N);//Initial field setting
	//datin();//Initial field input

//**** Simulation start ******************************
//Nstep = 200;
iout = -1;
start: ;

	if((((int)(time1) % Nstep)==0)) {datsave(ph, ND, N);}//Save the field every fixed repetition count
	if((((int)(time1) % Nstep)==0)) {datsave_paraview(ph, ND, N);}//Save the field every fixed repetition count
	//if(time1==200.) {datsave();}//save the field for a specific time

//**** Examine n00h[i][j] and m00h[n00][i][j] in each differential block *********************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ip=i+1; im=i-1;  jp=j+1;  jm=j-1;
			if(i==ndm){ip=0;}  if(i==0){im=ndm;}//periodic boundary conditions
			if(j==ndm){jp=0;}  if(j==0){jm=ndm;}

//--- Number of orientations where p is not 0 at location (i,j) and its surroundings (iÅ}1,jÅ}1)---
			n00=0;
			for(ii=1;ii<=nm;ii++){
				//if( (ph[ii][i][j]>0.0)||
				//	( (ph[ii][i][j]==0.0)&&(ph[ii][ip][j]>0.0)||
				//	  (ph[ii][im][j]>0.0)||
				//	  (ph[ii][i][jp]>0.0)||
				//	  (ph[ii][i][jm]>0.0)   ) ){
				//	  	n00++; m00h[n00][i][j]=ii;
				//		//printf("%d  ", n00);
				//}
				if( (ph[ii*ND*ND+i*ND+j]>0.0)||
					( (ph[ii*ND*ND+i*ND+j]==0.0)&&(ph[ii*ND*ND+ip*ND+j]>0.0)||
					  (ph[ii*ND*ND+im*ND+j]>0.0)||
					  (ph[ii*ND*ND+i*ND+jp]>0.0)||
					  (ph[ii*ND*ND+i*ND+jm]>0.0) ) ){
					  	n00++; m00h[n00*ND*ND+i*ND+j]=ii;
						//printf("%d  ", n00);
				}
			}
			//n00h[i][j]=n00;
			n00h[i*ND+j]=n00;
//--------------------------------------------------------------------------
		}
	}


//***** Evolution equation calculation **********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ip=i+1; im=i-1;  jp=j+1;  jm=j-1;
			if(i==ndm){ip=0;}  if(i==0){im=ndm;}//periodic boundary conditions
			if(j==ndm){jp=0;}  if(j==0){jm=ndm;}

			//for(n1=1; n1<=n00h[i][j]; n1++){
			for(n1=1; n1<=n00h[i*ND+j]; n1++){
				//ii=m00h[n1][i][j];  pddtt=0.0;
				ii=m00h[n1*ND*ND+i*ND+j];  pddtt=0.0;
				//for(n2=1; n2<=n00h[i][j]; n2++){
				for(n2=1; n2<=n00h[i*ND+j]; n2++){
					//jj=m00h[n2][i][j];  sum1=0.0;
					jj=m00h[n2*ND*ND+i*ND+j];  sum1=0.0;
					//for(n3=1; n3<=n00h[i][j]; n3++){
					for(n3=1; n3<=n00h[i*ND+j]; n3++){
						//kk=m00h[n3][i][j];
						kk=m00h[n3*ND*ND+i*ND+j];
						//sum1+=0.5*(aij[ii][kk]-aij[jj][kk])*(ph[kk][ip][j]+ph[kk][im][j]
						//									+ph[kk][i][jp]+ph[kk][i][jm]-4.0*ph[kk][i][j])
						//		 +(wij[ii][kk]-wij[jj][kk])*ph[kk][i][j];//[part of formula (4.31)]
						sum1+=0.5*(aij[ii*ND+kk]-aij[jj*ND+kk])*(ph[kk*ND*ND+ip*ND+j]+ph[kk*ND*ND+im*ND+j]
																+ph[kk*ND*ND+i*ND+jp]+ph[kk*ND*ND+i*ND+jm]-4.0*ph[kk*ND*ND+i*ND+j])
								 +(wij[ii*ND+kk]-wij[jj*ND+kk])*(ph[kk*ND*ND+i*ND+j]);	//[part of formula (4.31)]
					}
					//pddtt+=-2.0*tij[ii][jj]/double(n00h[i][j])
					//	*(sum1-8.0/PI*eij[ii][jj]*sqrt(ph[ii][i][j]*ph[jj][i][j]));
					pddtt+=-2.0*tij[ii*ND+jj]/double(n00h[i*ND+j])
				  *(sum1-8.0/PI*eij[ii*ND+jj]*sqrt(ph[ii*ND*ND+i*ND+j]*ph[jj*ND*ND+i*ND+j]));
					//Phase-field evolution equation [equation (4.31)]
				}
				//ph2[ii][i][j]=ph[ii][i][j]+pddtt*delt;		//Phase field time evolution (explicit method)
				//if(ph2[ii][i][j]>=1.0){ph2[ii][i][j]=1.0;}//Phase field domain correction
				//if(ph2[ii][i][j]<=0.0){ph2[ii][i][j]=0.0;}
				ph2[ii*ND*ND+i*ND+j]=ph[ii*ND*ND+i*ND+j]+pddtt*delt;		//Phase field time evolution (explicit method)
				if(ph2[ii*ND*ND+i*ND+j]>=1.0){ph2[ii*ND*ND+i*ND+j]=1.0;}//Phase field domain correction
				if(ph2[ii*ND*ND+i*ND+j]<=0.0){ph2[ii*ND*ND+i*ND+j]=0.0;}
			}
		}//j
	}//i

	for(k=1;k<=nm;k++){
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				//ph[k][i][j]=ph2[k][i][j];//move auxiliary array to main array
				ph[k*ND*ND+i*ND+j]=ph2[k*ND*ND+i*ND+j];	//move auxiliary array to main array
			}
		}
	}

//*** Phase field normalization correction ***********************************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//sum1=0.0; for(k=1;k<=nm;k++){ sum1+=ph[k][i][j]; }
			//for(k=1;k<=nm;k++){ ph[k][i][j]=ph[k][i][j]/sum1; }
			sum1=0.0; for(k=1;k<=nm;k++){ sum1+=ph[k*ND*ND+i*ND+j]; }
			for(k=1;k<=nm;k++){ ph[k*ND*ND+i*ND+j]=ph[k*ND*ND+i*ND+j]/sum1; }
		}
	}

//*********************************************************************
	//if(keypress()){return 0;}//Waiting for key
	time1=time1+1.;  if(time1<time1max){goto start;}//Determining if the maximum count has been reached

	end:;
  return 0;
}

//************ Initial field (phase field) setting subroutine *************
void ini000(double *ph, int ND, int N)
{
	int i, j, k, l, it;		//integer
	int ii, jj, kk;			//integer
	int ip, im, jp, jm;		//integer
	int x1, y1, x1h[10], y1h[10];//Coordinates of initial nuclei
	double sum1, t, r0, phi, r;
 	srand(3.0); // random number initialization
 	//srand(time(NULL)); // random number initialization
	int nd=ND, ndm=ND-1, nm=N-1, nmm=N-2;

	//x1h[1]=0.2*nd;   y1h[1]=0.2*nd;	//Coordinate setting of initial nucleus 1
	//x1h[2]=0.75*nd;  y1h[2]=0.4*nd;	//Coordinate setting of initial nucleus 2
	//x1h[3]=0.5*nd;   y1h[3]=0.75*nd;	//Coordinate setting of initial nucleus 3

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//for(ii=1;ii<=nm-1;ii++){ph[ii][i][j]=0.0;}
			//ph[nm][i][j]=1.0;//Initialize nmth phase field to 1
			for(ii=1;ii<=nm-1;ii++){ph[ii*ND*ND+i*ND+j]=0.0;}
			ph[nm*ND*ND+i*ND+j]=1.0;//Initialize nmth phase field to 1
		}
	}

	r0=5.0;
	for(ii=1;ii<=nm-1;ii++){
		//x1=x1h[ii]; y1=y1h[ii];
		x1=nd*DRND(1); y1=nd*DRND(1);//Initial nucleus position
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				r=sqrt( (double(i-x1))*(double(i-x1))+(double(j-y1))*(double(j-y1)) );
				//if(r<=r0){ ph[ii][i][j]=1.0;  ph[nm][i][j]=0.0; } //Set phase field for initial nuclear position
				if(r<=r0){
					ph[ii*ND*ND+i*ND+j]=1.0;
					ph[nm*ND*ND+i*ND+j]=0.0;
				} //Set phase field for initial nuclear position
			}
		}
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//sum1=0.0; for(k=1;k<=nm;k++){ sum1+=ph[k][i][j]; }
			//for(k=1;k<=nm;k++){ ph[k][i][j]=ph[k][i][j]/sum1; }//Phase field normalization correction
			sum1=0.0; for(k=1;k<=nm;k++){ sum1+=ph[k*ND*ND+i*ND+j]; }
			for(k=1;k<=nm;k++){ ph[k*ND*ND+i*ND+j]=ph[k*ND*ND+i*ND+j]/sum1; }//Phase field normalization correction
		}
	}

}

//************ data save subroutine *******************************
void datsave(double *ph, int ND, int N)
{
	FILE		*stream;		//Stream pointer setting
	int 		i, j, k;		//integer
	int nd=ND, ndm=ND-1, nm=N-1, nmm=N-2;

	stream = fopen("test.dat", "a");	//Open the file to write to in append mode
	fprintf(stream, "%e  \n", time1);	//Saving calculation counts
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=nm;k++){
				//fprintf(stream, "%e   ", ph[k][i][j]);//Save phase field
				fprintf(stream, "%e   ", ph[k*ND*ND+i*ND+j]);//Save phase field
			}
		}
	}
	fprintf(stream, "\n");	//writing a newline
	fclose(stream);					//close file
}

void datsave_paraview(double *ph, int ND, int N)
{
	FILE	*fp;
	char	fName[256];
	int 	i, j, k;
	int nd=ND, ndm=ND-1, nm=N-1, nmm=N-2;
	double *pht  = (double *)malloc(sizeof(double)*( ND*ND + ND ));
	
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			pht[i*ND+j]=0.0;
		}
	}
	
	iout = iout + 1;
	printf("paraview output no.%06d \n",iout);
	for(k=1;k<=nm;k++){
		sprintf(fName,"mpf1_N%03d_result%06d.vtk",k,iout);
		fp = fopen(fName, "w");
		fprintf(fp,"# vtk DataFile Version 3.0 \n");
		fprintf(fp,"output.vtk \n");
		fprintf(fp,"ASCII \n");
		fprintf(fp,"DATASET STRUCTURED_POINTS \n");
		fprintf(fp,"DIMENSIONS %5d %5d %5d \n",(ndm+1),(ndm+1),1);
		fprintf(fp,"ORIGIN 0.0 0.0 0.0 \n");
		fprintf(fp,"ASPECT_RATIO 1 1 1 \n");
		fprintf(fp,"POINT_DATA %16d \n",((ndm+1)*(ndm+1)*1));
		fprintf(fp,"SCALARS phase_field float \n");
		fprintf(fp,"LOOKUP_TABLE default \n");
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				//fprintf(stream, "%e   ", ph[k][i][j]);//Save phase field
				//fprintf(fp,"%10.6f\n", ph[k][i][j]);
				fprintf(fp,"%10.6f\n", ph[k*ND*ND+i*ND+j]);
				pht[i*ND+j]+=ph[k*ND*ND+i*ND+j]*float(k);
			}
		}
		fclose(fp);
	}
	
	for(k=0;k<=0;k++){
		sprintf(fName,"mpf1_N%03d_result%06d.vtk",k,iout);
		fp = fopen(fName, "w");
		fprintf(fp,"# vtk DataFile Version 3.0 \n");
		fprintf(fp,"output.vtk \n");
		fprintf(fp,"ASCII \n");
		fprintf(fp,"DATASET STRUCTURED_POINTS \n");
		fprintf(fp,"DIMENSIONS %5d %5d %5d \n",(ndm+1),(ndm+1),1);
		fprintf(fp,"ORIGIN 0.0 0.0 0.0 \n");
		fprintf(fp,"ASPECT_RATIO 1 1 1 \n");
		fprintf(fp,"POINT_DATA %16d \n",((ndm+1)*(ndm+1)*1));
		fprintf(fp,"SCALARS phase_field float \n");
		fprintf(fp,"LOOKUP_TABLE default \n");
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				fprintf(fp,"%10.6f\n", pht[i*ND+j]);
			}
		}
		fclose(fp);
	}
	free(pht);
}
//*********** data entry subroutine **************************
void datin(double *ph, int ND, int N)
{
	FILE		*datin0;//Stream pointer setting
	int 		i, j, k;//integer
	int nd=ND, ndm=ND-1, nm=N-1, nmm=N-2;

	datin0 = fopen("test.dat", "r");//open source file
	fscanf(datin0, "%lf", &time1);	//Read calculation count
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=nm;k++){
				//fscanf(datin0, "%lf  ", &ph[k][i][j]);//Load phase field
				fscanf(datin0, "%lf  ", &ph[k*ND*ND+i*ND+j]);//Load phase field
			}
		}
	}
	fclose(datin0);//close file

}

