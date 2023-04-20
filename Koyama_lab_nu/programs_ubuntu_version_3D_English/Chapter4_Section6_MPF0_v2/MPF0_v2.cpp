//Program for polycrystalline grain structure formation
//Grain number is 1-nm
//grain number nm liquid phase

#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

#define DRND(x) ((double)(x)/RAND_MAX*rand())//Random number setting

//#define ND 100			//Number of divisions per side of calculation domain in difference calculation
//#define N 5				//Number of crystal orientations to consider + 1

	//int nd=ND, ndm=ND-1;	//Define the number of difference divisions
							// (number of difference blocks) on one side of the computational domain, ND-1
	//int nm=N-1, nmm=N-2;	//Define the number of crystal orientations to consider,
							// N-2 (number of crystal orientations to consider - 1)
	double PI=3.141592;		//pi
	double RR=8.3145;		//gas constant
	double time1;			//calculation count number
	//double ph[N][ND][ND][ND], ph2[N][ND][ND][ND];	//phase field, phase field auxiliary array
	//double aij[N][N];//gradient energy factor
	//double wij[N][N];//the coefficient of the penalty term
	//double tij[N][N];//Grain boundary mobility
	//double eij[N][N];//Driving force of grain boundary migration
	//int m00h[N][ND][ND][ND];//number of directions where p is not 0 at location (i,j) and its surroundings (iÅ}1,jÅ}1)
	//int n00h[ND][ND][ND];//Number of orientations where p is not 0 at location (i,j) and its surroundings (iÅ}1,jÅ}1)
	int Nstep, iout;

	void ini000(double *ph, int ND, int N);	//Initial field setup subroutine
	void datin(double *ph, int ND, int N);	//Subroutine for initial field reading
	void datsave(double *ph, int ND, int N);	//data save subroutine
	void datsave_paraview(double *ph, int ND, int N);	//data save subroutine

//******* main program ******************************************
int main(void)
{
	int ND, N;
	int nd, ndm, nm, nmm;

	int i, j, k, l, ii, jj, kk, ll, it;	//integer
	int ip, im, jp, jm, kp, km;			//integer
	int n1, n2, n3;						//integer
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

	double dtemp;
	int    readff;

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
	dtemp   = data[4];
	L       = data[5];
	vm0     = data[6];
	gamma0  = data[7];
	delta   = data[8];
	K0      = data[9];
	W0      = data[10];
	amobi   = data[11];
	M0      = data[12];
	E0      = data[13];
	time1max= int(data[14]);
	Nstep   = int(data[15]);
	readff  = int(data[16]);
	printf("---------------------------------\n");
	//
	nd=ND;
	ndm=ND-1;
	nm=N-1;
	nmm=N-2;
	//
	double *ph   = (double *)malloc(sizeof(double)*( N*ND*ND*ND ));	//phase field
	double *ph2  = (double *)malloc(sizeof(double)*( N*ND*ND*ND ));	//Phase field auxiliary array
	double *aij  = (double *)malloc(sizeof(double)*( ND*ND ));	//gradient energy factor
	double *wij  = (double *)malloc(sizeof(double)*( ND*ND ));	//the coefficient of the penalty term
	double *tij  = (double *)malloc(sizeof(double)*( ND*ND ));	//Grain boundary mobility
	double *eij  = (double *)malloc(sizeof(double)*( ND*ND ));	//Driving force of grain boundary migration
	double *m00h = (double *)malloc(sizeof(double)*( N*ND*ND*ND ));	//number of directions where p is not 0 at location (i,j) and its surroundings (iÅ}1,jÅ}1)
	int    *n00h =       (int *)malloc(sizeof(int)*( ND*ND*ND ));	//Number of orientations where p is not 0 at location (i,j) and its surroundings (iÅ}1,jÅ}1)
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

	time1=0.;					//Initial value of calculation count
	//time1max=10000001.0;		//Maximum calculation count

//*** Equation (4.32) - set array (K,W,M,E) in equation (4.35) *****************************
	for(i=1;i<=nm;i++){
		for(j=1;j<=nm;j++){
			//wij[i][j]=W1;  aij[i][j]=K1;  tij[i][j]=0.0;  eij[i][j]=0.0;
			//if( (i==nm)||(j==nm) ){eij[i][j]=E1; tij[i][j]=M1;}
			//if(i>j){eij[i][j]=-eij[i][j];}
			//if(i==j){	wij[i][j]=0.0;  aij[i][j]=0.0;  tij[i][j]=0.0; eij[i][j]=0.0;}
			wij[i*ND+j]=W1;  aij[i*ND+j]=K1;  tij[i*ND+j]=0.0;  eij[i*ND+j]=0.0;
			if( (i==nm)||(j==nm) ){eij[i*ND+j]=E1; tij[i*ND+j]=M1;}
			if(i>j){eij[i*ND+j]=-eij[i*ND+j];}
			if(i==j){wij[i*ND+j]=0.0;  aij[i*ND+j]=0.0;  tij[i*ND+j]=0.0; eij[i*ND+j]=0.0;}
		}
	}

//*** Initial Field Settings and Drawing Window Display *****************************************
	for(kk=1;kk<=nm;kk++){
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				for(k=0;k<=ndm;k++){
					//ph[kk][i][j][k]=0.0;  ph2[kk][i][j][k]=0.0;//Phase field and auxiliary array initialization
					 ph[kk*ND*ND*ND+i*ND*ND+j*ND+k]=0.0;	//Phase field initialization
					ph2[kk*ND*ND*ND+i*ND*ND+j*ND+k]=0.0;	//Auxiliary array initialization
				}
			}
		}
	}
	
	if(readff==0){
		printf("make initial concentration (random) \n");
		ini000(ph, ND, N);	//Initial field setting
	} else {
		printf("read data.dat file \n");
		datin(ph, ND, N);
	}
//**** Simulation start ******************************
//Nstep = 20;
iout = -1;
start: ;

	if((((int)(time1) % Nstep)==0)) {iout = iout + 1;}
	if((((int)(time1) % Nstep)==0)) {datsave(ph, ND, N);}	//Save the field every fixed repetition count
	if((((int)(time1) % Nstep)==0)) {datsave_paraview(ph, ND, N);}	//Save the field every fixed repetition count
	//if(time1==200.) {datsave();}					//save the field for a specific time

//**** Examine n00h[i][j] and m00h[n00][i][j] in each differential block *********************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				
				//periodic boundary conditions
				ip=i+1; im=i-1;  jp=j+1;  jm=j-1;  kp=k+1;  km=k-1;
				if(i==ndm){ip=0;}  if(i==0){im=ndm;}
				if(j==ndm){jp=0;}  if(j==0){jm=ndm;}
				if(k==ndm){kp=0;}  if(k==0){km=ndm;}

//--- Number of orientations where p is not 0 at location (i,j) and its surroundings (iÅ}1,jÅ}1)---
				n00=0;
				for(ii=1;ii<=nm;ii++){
					//if( (ph[ii][i][j][k]>0.0)||
					//  ( (ph[ii][i][j][k]==0.0)&&
					//		(ph[ii][ip][j][k]>0.0)||ph[ii][im][j][k]>0.0)||
					//  	(ph[ii][i][jp][k]>0.0)||ph[ii][i][jm][k]>0.0)||
					//		(ph[ii][i][j][kp]>0.0)||ph[ii][i][j][km]>0.0) ) ){
					//  		n00++; m00h[n00][i][j][k]=ii;
					//}
					if( (ph[ii*ND*ND*ND+i*ND*ND+j*ND+k]>0.0)||
						( (ph[ii*ND*ND*ND+i*ND*ND+j*ND+k]==0.0)&&(
							(ph[ii*ND*ND*ND+ip*ND*ND+j*ND+k]>0.0) || (ph[ii*ND*ND*ND+im*ND*ND+j*ND+k]>0.0) ||
							(ph[ii*ND*ND*ND+i*ND*ND+jp*ND+k]>0.0) || (ph[ii*ND*ND*ND+i*ND*ND+jm*ND+k]>0.0) ||
							(ph[ii*ND*ND*ND+i*ND*ND+j*ND+kp]>0.0) || (ph[ii*ND*ND*ND+i*ND*ND+j*ND+km]>0.0) )
						) ){n00++; m00h[n00*ND*ND*ND+i*ND*ND+j*ND+k]=ii;}
				}
				//n00h[i][j][k]=n00;
				n00h[i*ND*ND+j*ND+k]=n00;
//--------------------------------------------------------------------------
			}
		}
	}

//***** Evolution equation calculation **********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){

				//periodic boundary conditions
				ip=i+1; im=i-1;  jp=j+1;  jm=j-1;  kp=k+1;  km=k-1;
				if(i==ndm){ip=0;}  if(i==0){im=ndm;}
				if(j==ndm){jp=0;}  if(j==0){jm=ndm;}
				if(k==ndm){kp=0;}  if(k==0){km=ndm;}

				//for(n1=1; n1<=n00h[i][j][k]; n1++){
				for(n1=1; n1<=n00h[i*ND*ND+j*ND+k]; n1++){
					//ii=m00h[n1][i][j][k];  pddtt=0.0;
					ii=m00h[n1*ND*ND*ND+i*ND*ND+j*ND+k];  pddtt=0.0;
					//for(n2=1; n2<=n00h[i][j][k]; n2++){
					for(n2=1; n2<=n00h[i*ND*ND+j*ND+k]; n2++){
						//jj=m00h[n2][i][j][k];  sum1=0.0;
						jj=m00h[n2*ND*ND*ND+i*ND*ND+j*ND+k];  sum1=0.0;
						//for(n3=1; n3<=n00h[i][j][k]; n3++){
						for(n3=1; n3<=n00h[i*ND*ND+j*ND+k]; n3++){
							//kk=m00h[n3][i][j][k];
							kk=m00h[n3*ND*ND*ND+i*ND*ND+j*ND+k];
							//sum1+=0.5*(aij[ii][kk]-aij[jj][kk])*(ph[kk][ip][j][k]+ph[kk][im][j][k]
							//									  +ph[kk][i][jp][k]+ph[kk][i][jm][k]
							//									  +ph[kk][i][j][kp]+ph[kk][i][j][km]
							//									  -6.0*ph[kk][i][j][k])
							//		 +(wij[ii][kk]-wij[jj][kk])*ph[kk][i][j][k];	//[part of formula (4.31)]
							sum1+=0.5*(aij[ii*ND+kk]-aij[jj*ND+kk])*(ph[kk*ND*ND*ND+ip*ND*ND+j*ND+k]+ph[kk*ND*ND*ND+im*ND*ND+j*ND+k]
																	+ph[kk*ND*ND*ND+i*ND*ND+jp*ND+k]+ph[kk*ND*ND*ND+i*ND*ND+jm*ND+k]
																	+ph[kk*ND*ND*ND+i*ND*ND+j*ND+kp]+ph[kk*ND*ND*ND+i*ND*ND+j*ND+km]
																	-6.0*ph[kk*ND*ND*ND+i*ND*ND+j*ND+k])
									 +(wij[ii*ND+kk]-wij[jj*ND+kk])*(ph[kk*ND*ND*ND+i*ND*ND+j*ND+k]);	//[part of formula (4.31)]
						}
						//pddtt+=-2.0*tij[ii][jj]/double(n00h[i][j][k])
						//	*(sum1-8.0/PI*eij[ii][jj]*sqrt(ph[ii][i][j][k]*ph[jj][i][j][k]));
						pddtt+=-2.0*tij[ii*ND+jj]/double(n00h[i*ND*ND+j*ND+k])
							  *(sum1-8.0/PI*eij[ii*ND+jj]*sqrt(ph[ii*ND*ND*ND+i*ND*ND+j*ND+k]*ph[jj*ND*ND*ND+i*ND*ND+j*ND+k]));
						//Phase-field evolution equation [equation (4.31)]
					}
					//ph2[ii][i][j][k]=ph[ii][i][j][k]+pddtt*delt;		//Phase field time evolution (explicit method)
					//if(ph2[ii][i][j][k]>=1.0){ph2[ii][i][j][k]=1.0;}//Phase field domain correction
					//if(ph2[ii][i][j][k]<=0.0){ph2[ii][i][j][k]=0.0;}
					ph2[ii*ND*ND*ND+i*ND*ND+j*ND+k]=ph[ii*ND*ND*ND+i*ND*ND+j*ND+k]+pddtt*delt;		//Phase field time evolution (explicit method)
					if(ph2[ii*ND*ND*ND+i*ND*ND+j*ND+k]>=1.0){ph2[ii*ND*ND*ND+i*ND*ND+j*ND+k]=1.0;}//Phase field domain correction
					if(ph2[ii*ND*ND*ND+i*ND*ND+j*ND+k]<=0.0){ph2[ii*ND*ND*ND+i*ND*ND+j*ND+k]=0.0;}
				}
			}//k
		}//j
	}//i

	for(kk=1;kk<=nm;kk++){
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				for(k=0;k<=ndm;k++){
					//ph[kk][i][j][k]=ph2[kk][i][j][k];	//move auxiliary array to main array
					ph[kk*ND*ND*ND+i*ND*ND+j*ND+k]=ph2[kk*ND*ND*ND+i*ND*ND+j*ND+k];	//move auxiliary array to main array
				}
			}
		}
	}

//*** Phase field normalization correction ***********************************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//sum1=0.0; for(kk=1;kk<=nm;kk++){ sum1+=ph[k][i][j][k]; }
				//for(kk=1;kk<=nm;kk++){ ph[k][i][j][k]=ph[k][i][j][k]/sum1; }
				sum1=0.0;
				for(kk=1;kk<=nm;kk++){ sum1+=ph[kk*ND*ND*ND+i*ND*ND+j*ND+k]; }
				for(kk=1;kk<=nm;kk++){ ph[kk*ND*ND*ND+i*ND*ND+j*ND+k]=ph[kk*ND*ND*ND+i*ND*ND+j*ND+k]/sum1; }
			}
		}
	}

//*********************************************************************
	//if(keypress()){return 0;}//Waiting for key
	time1=time1+1.;
	temp += dtemp * delt;
	if(time1<time1max){goto start;}//Determining if the maximum count has been reached

end:;
	printf("Finished \n");
	std::exit(0);
	//return 0;
}

//************ Initial field (phase field) setting subroutine *************
void ini000(double *ph, int ND, int N)
{
	int i, j, k, l, it;		//integer
	int ii, jj, kk;			//integer
	int ip, im, jp, jm, kp, km;	//integer
	int x1, y1, z1, x1h[10], y1h[10], z1h[10];//Coordinates of initial nuclei
	double sum1, t, r0, phi, r;
	//srand(time(NULL)); // random number initialization
	int nd=ND, ndm=ND-1, nm=N-1, nmm=N-2;

	//x1h[1]=0.2*nd;   y1h[1]=0.2*nd;   z1h[1]=0.2*nd;	//Coordinate setting of initial nucleus 1
	//x1h[2]=0.75*nd;  y1h[2]=0.4*nd;   z1h[2]=0.4*nd;	//Coordinate setting of initial nucleus 2
	//x1h[3]=0.5*nd;   y1h[3]=0.75*nd;  z1h[3]=0.75*nd;	//Coordinate setting of initial nucleus 3

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//for(ii=1;ii<=nm-1;ii++){ph[ii][i][j][k]=0.0;}
				//ph[nm][i][j][k]=1.0;//Initialize nmth phase field to 1
				for(ii=1;ii<=nm-1;ii++){ph[ii*ND*ND*ND+i*ND*ND+j*ND+k]=0.0;}
				ph[nm*ND*ND*ND+i*ND*ND+j*ND+k]=1.0;//Initialize nmth phase field to 1
			}
		}
	}

	r0=5.0;
	for(ii=1;ii<=nm-1;ii++){
		//x1=x1h[ii]; y1=y1h[ii]; z1=z1h[ii];
		x1=nd*DRND(1); y1=nd*DRND(1); z1=nd*DRND(1);
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				for(k=0;k<=ndm;k++){
					r=sqrt( (double(i-x1))*(double(i-x1))+(double(j-y1))*(double(j-y1))+(double(k-z1))*(double(k-z1)) );
					//if(r<=r0){ ph[ii][i][j][k]=1.0;  ph[nm][i][j][k]=0.0; } //Set phase field for initial nuclear position
					if(r<=r0){
						ph[ii*ND*ND*ND+i*ND*ND+j*ND+k]=1.0;
						ph[nm*ND*ND*ND+i*ND*ND+j*ND+k]=0.0;
					} //Set phase field for initial nuclear position
				}
			}
		}
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//sum1=0.0; for(k=1;k<=nm;k++){ sum1+=ph[kk][i][j][k]; }
				//for(k=1;k<=nm;k++){ ph[kk][i][j][k]=ph[kk][i][j][k]/sum1; }//Phase field normalization correction
				sum1=0.0;
				for(kk=1;kk<=nm;kk++){ sum1+=ph[kk*ND*ND*ND+i*ND*ND+j*ND+k]; }
				for(kk=1;kk<=nm;kk++){ ph[kk*ND*ND*ND+i*ND*ND+j*ND+k]=ph[kk*ND*ND*ND+i*ND*ND+j*ND+k]/sum1; }//Phase field normalization correction
			}
		}
	}

}

//************ data save subroutine *******************************
void datsave(double *ph, int ND, int N)
{
	FILE	*stream;		//Stream pointer setting
	char	fName[256];
	int 	i, j, k, kk;		//integer
	int nd=ND, ndm=ND-1, nm=N-1, nmm=N-2;
	int ndxm=ndm, ndym=ndm, ndzm=ndm;

	sprintf(fName,"data_%06d.dat",iout);
	//stream = fopen("test.dat", "a");	//Open the file to write to in append mode
	stream = fopen(fName, "w");
	fprintf(stream, "%d %d %d \n", ndm, ndm, ndm);
	fprintf(stream, "%e  \n", time1);	//Saving calculation counts
	for(i=0;i<=ndxm;i++){
		for(j=0;j<=ndym;j++){
			for(k=0;k<=ndzm;k++){
				for(kk=0;kk<=nm;kk++){
					//fprintf(stream, "%e   ", ph[kk][i][j][k]);//Save phase field
					fprintf(stream, "%e   ", ph[kk*ND*ND*ND+i*ND*ND+j*ND+k]);//Save phase field
				}
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
	int 	i, j, k, kk;
	int nd=ND, ndm=ND-1, nm=N-1, nmm=N-2;
	int ndxm=ndm, ndym=ndm, ndzm=ndm;
	double *pht  = (double *)malloc(sizeof(double)*( ND*ND*ND ));
	
	for(i=0;i<=ndxm;i++){
		for(j=0;j<=ndym;j++){
			for(k=0;k<=ndzm;k++){
				pht[i*ND*ND+j*ND+k]=0.0;
			}
		}
	}
	
	printf("paraview output no.%06d \n",iout);
	for(kk=1;kk<=nm;kk++){
		sprintf(fName,"mpf0_N%03d_result%06d.vtk",kk,iout);
		fp = fopen(fName, "w");
		fprintf(fp,"# vtk DataFile Version 3.0 \n");
		fprintf(fp,"output.vtk \n");
		fprintf(fp,"ASCII \n");
		fprintf(fp,"DATASET STRUCTURED_POINTS \n");
		fprintf(fp,"DIMENSIONS %5d %5d %5d \n",(ndxm+1),(ndym+1),(ndzm+1));
		fprintf(fp,"ORIGIN 0.0 0.0 0.0 \n");
		fprintf(fp,"ASPECT_RATIO %f %f %f \n",float(ndxm/ndxm),float(ndym/ndxm),float(ndzm/ndxm));
		fprintf(fp,"POINT_DATA %16d \n",((ndxm+1)*(ndym+1)*(ndzm+1)));
		fprintf(fp,"SCALARS phase_field float \n");
		fprintf(fp,"LOOKUP_TABLE default \n");
		for(k=0;k<=ndzm;k++){
			for(j=0;j<=ndym;j++){
				for(i=0;i<=ndxm;i++){
					//fprintf(stream, "%e   ", ph[kk][i][j][k]);//Save phase field
					//fprintf(fp,"%10.6f\n", ph[kk][i][j][k]);
					fprintf(fp,"%10.6f\n", ph[kk*ND*ND*ND+i*ND*ND+j*ND+k]);
					pht[i*ND*ND+j*ND+k]+=ph[kk*ND*ND*ND+i*ND*ND+j*ND+k]*float(kk);
				}
			}
		}
		fclose(fp);
	}
	
	for(kk=0;kk<=0;kk++){
		sprintf(fName,"mpf0_N%03d_result%06d.vtk",kk,iout);
		fp = fopen(fName, "w");
		fprintf(fp,"# vtk DataFile Version 3.0 \n");
		fprintf(fp,"output.vtk \n");
		fprintf(fp,"ASCII \n");
		fprintf(fp,"DATASET STRUCTURED_POINTS \n");
		fprintf(fp,"DIMENSIONS %5d %5d %5d \n",(ndxm+1),(ndym+1),(ndzm+1));
		fprintf(fp,"ORIGIN 0.0 0.0 0.0 \n");
		fprintf(fp,"ASPECT_RATIO %f %f %f \n",float(ndxm/ndxm),float(ndym/ndxm),float(ndzm/ndxm));
		fprintf(fp,"POINT_DATA %16d \n",((ndxm+1)*(ndym+1)*(ndzm+1)));
		fprintf(fp,"SCALARS phase_field float \n");
		fprintf(fp,"LOOKUP_TABLE default \n");
		for(k=0;k<=ndzm;k++){
			for(j=0;j<=ndym;j++){
				for(i=0;i<=ndxm;i++){
					fprintf(fp,"%10.6f\n", pht[i*ND*ND+j*ND+k]);
				}
			}
		}
		fclose(fp);
	}
	free(pht);
}
//*********** data entry subroutine **************************
void datin(double *ph, int ND, int N)
{
	FILE	*datin0;//Stream pointer setting
	int 	i, j, k, kk;//integer
	int nd=ND, ndm=ND-1, nm=N-1, nmm=N-2;
	int ndxm=ndm, ndym=ndm, ndzm=ndm;

	datin0 = fopen("data.dat", "r");//Open the source file

	fscanf(datin0, "%d %d %d", &rndxm, &rndym, &rndzm);
	if (ndm != rndxm){
		printf("data size is mismatch \n");
		printf("Please, change ND= %d in parameters.txt \n", rndxm);
	}
	
	fscanf(datin0, "%lf", &time1);
	
	for(i=0;i<=ndxm;i++){
		for(j=0;j<=ndym;j++){
			for(k=0;k<=ndzm;k++){
				for(kk=0;kk<=nm;kk++){
					fscanf(datin0, "%lf  ", &ph[kk*ND*ND*ND+i*ND*ND+j*ND+k]);//Load phase field
				}
			}
		}
	}
	printf("time=  %f  \n", time1);
	fclose(datin0);//close file

}