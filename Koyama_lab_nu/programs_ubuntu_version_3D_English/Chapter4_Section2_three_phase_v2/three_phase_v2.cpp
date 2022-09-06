#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

#define DRND(x) ((double)(x)/RAND_MAX*rand())  //Random number setting

//#define ND 100					//Number of divisions per side of calculation domain in difference calculation

	//int nd=ND;					//Number of difference divisions on one side of the computational domain (number of difference blocks)
	//int ndm=ND-1;				//Define ND-1
	double PI=3.141592654;		//Pi
	double rr=8.3145;			//gas constant
	double temp;				//absolute temperature
	double time1;				//Number of calculation counts (proportional to time)
	double c2a, c3a;			//Average composition (1: A component, 2: B component, 3: C component)
	//double c2h[ND][ND][ND, c3h[ND][ND][ND];//topical composition
	int Nstep, iout;

	void ini000(double *c2h, double *c3h, int ND);	//Initial concentration profile setting subroutine
	void datin(double *c2h, double *c3h, int ND);		//Subroutine for initial field reading
	void datsave(double *c2h, double *c3h, int ND);	//data save subroutine
	void datsave_paraview(double *c2h, double *c3h, int ND);//data save subroutine

//******* main program ******************************************
int main(void)
{
	int ND, nd, ndm;
	
	double c1, c2, c3;							//local concentration
	//double c2h2[ND][ND][ND], c3h2[ND][ND][ND];//auxiliary matrix for the local concentration field
	double c2k_chem, c2k_su;					//local potential
	//double c2k[ND][ND][ND];					//local potential
	double c3k_chem, c3k_su;					//local potential
	//double c3k[ND][ND][ND];					//local potential
	double dc2a, sumc2, dc3a, sumc3, sumc23;	//Variables used for concentration field balance calculation
	double dakd2, dakd3;						//Second derivative of the diffusion potential
	double c2ddtt, c3ddtt;						//Time variation of concentration field

	int   i, j, k;								//integer
	int   ip, im, jp, jm, kp, km;				//(i+1),(i-1),(j+1),(j-1),(k+1),(k-1)
	double al, b1, rtemp, delt;					//Length of one side of computational domain, length of one side of difference block, RT, time step
	double time1max;							//Maximum calculation count (calculation end count)
	double cmob22, cmob33, cmob23, cmob32;		//mobility

	double om_12, om_23, om_13;					//interaction parameter
	double om_12e, om_23e, om_13e;				//interaction parameter
	double kapa_c1, kapa_c2, kapa_c3;					//concentration gradient energy coefficient
	double kapa_c1c, kapa_c2c, kapa_c3c;					//concentration gradient energy coefficient
	
	int    readff;

//****** Setting calculation conditions and material constants ****************************************
	printf("---------------------------------\n");
	printf("read parameters from parameters.txt\n");
	FILE *fp;
	char name[40], comment[72];
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
	c2a     = data[1];
	c3a     = data[2];
	delt    = data[3];
	temp    = data[4];
	al      = data[5];	// [nm]
	cmob22  = data[6];
	cmob33  = data[7];
	cmob23  = data[8];
	om_12e  = data[9];	//L_AB
	om_13e  = data[10];	//L_AC
	om_23e  = data[11];	//L_BC
	kapa_c1c= data[12];
	kapa_c2c= data[13];
	kapa_c3c= data[14];
	time1max= int(data[15]);
	Nstep   = int(data[16]);
	readff  = int(data[17]);
	printf("---------------------------------\n");
	//
	nd=ND;					//Number of difference divisions on one side of the computational domain (number of difference blocks)
	ndm=ND-1;				//Define ND-1
	//
	//double ch[(nstep+1)*nstep]; // [(nstep+1)*nstep]=[nstep*nstep+nstep], a[i][j]= a[i*n+j] for a[][n]
	//a[z * ySize * xSize + y * xSize + x], a[i][j][k]=a[i*n*m + j*n + k]
	double *c2h  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//topical composition
	double *c3h  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//topical composition
	double *c2h2 = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//auxiliary matrix for the local concentration field
	double *c3h2 = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//auxiliary matrix for the local concentration field
	double *c2k  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//local potential
	double *c3k  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//local potential
	//
	//printf("C2B(0.333) =  ");	scanf(" %lf",&c2a);//Input average composition (B component) from standard input/output	//	c2a=1./3.;
	//printf("C3C(0.333) =  ");	scanf(" %lf",&c3a);//Input average composition (C component) from standard input/output	//	c3a=1./3.;
	//printf("delt(0.005)=  ");	scanf(" %lf",&delt);//ticking	//	delt=0.005;

	//temp=900.0;					//temperature (K)
	rtemp=rr*temp;				//RT
	//al=100.0*1.0E-09;			//Side length of calculation area (m)
	al=al*1.0E-09;			//Side length of calculation area (m)
	b1=al/(double)ND;			//Length of one side of difference block

	//cmob22=1.0;					//mobility
	//cmob33=1.0;					//mobility
	//cmob23=cmob32=-0.5;			//mobility
	cmob32=cmob23;

	//om_12=25000./rtemp; 		//Interaction parameters (in J/mol, non-dimensionalized at RT)
	//om_13=25000./rtemp;
	//om_23=25000./rtemp;
	om_12=om_12e/rtemp; 		//Interaction parameters (in J/mol, non-dimensionalized at RT)
	om_13=om_13e/rtemp;
	om_23=om_23e/rtemp;

	//kapa_c2=5.0e-15/b1/b1/rtemp;//gradient energy coefficient (in Jm^2/mol, non-dimensionalized by RT and b1^2)
	//kapa_c3=5.0e-15/b1/b1/rtemp;//gradient energy coefficient (in Jm^2/mol, non-dimensionalized by RT and b1^2)
	kapa_c1=kapa_c1c/b1/b1/rtemp;//gradient energy coefficient (in Jm^2/mol, non-dimensionalized by RT and b1^2)
	kapa_c2=kapa_c2c/b1/b1/rtemp;//gradient energy coefficient (in Jm^2/mol, non-dimensionalized by RT and b1^2)
	kapa_c3=kapa_c3c/b1/b1/rtemp;//gradient energy coefficient (in Jm^2/mol, non-dimensionalized by RT and b1^2)

	time1=0.0;					//Initial value of calculation count
	//time1max=1.0e05+1.0;		//Maximum calculation count

//*** Initial concentration field setting and drawing window display *****************************************
	if(readff==0){
		printf("make initial concentration (random) \n");
		ini000(c2h, c3h, ND);//Setting the initial concentration field
	} else {
		printf("read data.dat file \n");
		datin(c2h, c3h, ND);
	}
//**** Simulation start ******************************
//Nstep = 2000;
iout = -1;
start: ;

	//printf("Time: %f \n", time1);
	if((((int)(time1) % Nstep)==0)) {iout = iout + 1;}
	if((((int)(time1) % Nstep)==0)) {datsave(c2h, c3h, ND);} //Save the concentration field every fixed repetition count
	if((((int)(time1) % Nstep)==0)) {datsave_paraview(c2h, c3h, ND);} //Save the concentration field every fixed repetition count

//***** Potential field calculation ***********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				
				ip=i+1; im=i-1; jp=j+1; jm=j-1; kp=k+1; km=k-1;
				if(i==ndm){ip=0;} if(i==0){im=ndm;}	//periodic boundary conditions
				if(j==ndm){jp=0;} if(j==0){jm=ndm;}	//periodic boundary conditions
				if(k==ndm){kp=0;} if(k==0){km=ndm;}	//periodic boundary conditions

				//c2=c2h[i][j]; c3=c3h[i][j]; c1=1.0-c2-c3;//local concentration field
				c2=c2h[i*ND*ND+j*ND+k]; c3=c3h[i*ND*ND+j*ND+k]; c1=1.0-c2-c3;//local concentration field

				//Gsys = Gc + Ggrad
				//Gc = integral { (L_AB * c1 * c2) + (L_AC * c1 * c3) + (L_BC * c2 * c3) + ...
				//  + RT*( c1*log(c1) + c2*log(c2) + c3*log(c3) + ...) }  dr
				c2k_chem=om_12*(c1-c2)-om_13*c3+om_23*c3+(log(c2)-log(c1));//chemical diffusion potential
				//c2k_su=-2.*kapa_c2*(c2h[ip][j][k]+c2h[im][j][k]
      			//					 +c2h[i][jp][k]+c2h[i][jm][k]
				//					 +c2h[i][j][kp]+c2h[i][j][km]-6.0*c2)
				//			-kapa_c3*(c3h[ip][j][k]+c3h[im][j][k]
      			//					 +c3h[i][jp][k]+c3h[i][jm][k]
				//					 +c3h[i][j][kp]+c3h[i][j][km]-6.0*c3);//gradient potential

				//Ggrad = integral {0.5*kappa1*div(div c1) + 0.5*kappa2*div(div c2) + 0.5*kappa3*div(div c3) + ... } dr
				c2k_su=-2.0*(0.5*kapa_c1+0.5*kapa_c2)*(c2h[ip*ND*ND+j*ND+k]+c2h[im*ND*ND+j*ND+k]
													  +c2h[i*ND*ND+jp*ND+k]+c2h[i*ND*ND+jm*ND+k]
      												  +c2h[i*ND*ND+j*ND+kp]+c2h[i*ND*ND+j*ND+km]
													  -6.0*c2)
						   -(0.5*kapa_c1+0.5*kapa_c3)*(c3h[ip*ND*ND+j*ND+k]+c3h[im*ND*ND+j*ND+k]
													  +c3h[i*ND*ND+jp*ND+k]+c3h[i*ND*ND+jm*ND+k]
      												  +c3h[i*ND*ND+j*ND+kp]+c3h[i*ND*ND+j*ND+km]
													  -6.0*c3);	//gradient potential

				c3k_chem=om_13*(c1-c3)-om_12*c2+om_23*c2+(log(c3)-log(c1));///chemical diffusion potential
				//c3k_su=-2.*kapa_c3*(c3h[ip][j][k]+c3h[im][j][k]
      			//					 +c3h[i][jp][k]+c3h[i][jm][k]
				//					 +c3h[i][j][kp]+c3h[i][j][km]-6.0*c3)
				//			-kapa_c2*(c2h[ip][j][k]+c2h[im][j][k]
      			//					 +c2h[i][jp][k]+c2h[i][jm][k]
				//					 +c2h[i][j][kp]+c2h[i][j][km]-6.0*c2);//gradient potential
				
				c3k_su=-2.0*(0.5*kapa_c1+0.5*kapa_c3)*(c3h[ip*ND*ND+j*ND+k]+c3h[im*ND*ND+j*ND+k]
													  +c3h[i*ND*ND+jp*ND+k]+c3h[i*ND*ND+jm*ND+k]
      												  +c3h[i*ND*ND+j*ND+kp]+c3h[i*ND*ND+j*ND+km]
													  -6.0*c3)
						   -(0.5*kapa_c1+0.5*kapa_c2)*(c2h[ip*ND*ND+j*ND+k]+c2h[im*ND*ND+j*ND+k]
													  +c2h[i*ND*ND+jp*ND+k]+c2h[i*ND*ND+jm*ND+k]
      												  +c2h[i*ND*ND+j*ND+kp]+c2h[i*ND*ND+j*ND+km]
													  -6.0*c2);	//gradient potential
				
				//c2k[i][j][k]=c2k_chem+c2k_su;//Diffusion potential (equation (4.1))
				//c3k[i][j][k]=c3k_chem+c3k_su;
				c2k[i*ND*ND+j*ND+k]=c2k_chem+c2k_su;//Diffusion potential (equation (4.1))
				c3k[i*ND*ND+j*ND+k]=c3k_chem+c3k_su;
			}
		}
	}

//***** Evolution equation calculation **********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){

				ip=i+1; im=i-1; jp=j+1; jm=j-1; kp=k+1; km=k-1;
				if(i==ndm){ip=0;} if(i==0){im=ndm;}	//periodic boundary conditions
				if(j==ndm){jp=0;} if(j==0){jm=ndm;}	//periodic boundary conditions
				if(k==ndm){kp=0;} if(k==0){km=ndm;}	//periodic boundary conditions

				//dakd2=c2k[ip][j][k]+c2k[im][j][k]
      			//	   +c2k[i][jp][k]+c2k[i][jm][k]
				//	   +c2k[i][j][kp]+c2k[i][j][km]
				//	   -6.0*c2k[i][j][k]);
				//dakd3=c3k[ip][j][k]+c3k[im][j][k]
      			//	   +c3k[i][jp][k]+c3k[i][jm][k]
				//	   +c3k[i][j][kp]+c3k[i][j][km]
				//	   -6.0*c3k[i][j][k]);

				dakd2=c2k[ip*ND*ND+j*ND+k]+c2k[im*ND*ND+j*ND+k]
					 +c2k[i*ND*ND+jp*ND+k]+c2k[i*ND*ND+jm*ND+k]
      				 +c2k[i*ND*ND+j*ND+kp]+c2k[i*ND*ND+j*ND+km]
					 -6.0*c2k[i*ND*ND+j*ND+k];	//Second derivative of the diffusion potential
				dakd3=c3k[ip*ND*ND+j*ND+k]+c3k[im*ND*ND+j*ND+k]
					 +c3k[i*ND*ND+jp*ND+k]+c3k[i*ND*ND+jm*ND+k]
      				 +c3k[i*ND*ND+j*ND+kp]+c3k[i*ND*ND+j*ND+km]
					 -6.0*c3k[i*ND*ND+j*ND+k];	//Second derivative of the diffusion potential

				c2ddtt=cmob22*dakd2+cmob23*dakd3;//Diffusion equation (equation (4.2))
				c3ddtt=cmob32*dakd2+cmob33*dakd3;

				//c2h2[i][j][k]=c2h[i][j][k]+c2ddtt*delt;//Time evolution of the concentration field
				//c3h2[i][j][k]=c3h[i][j][k]+c3ddtt*delt;
				c2h2[i*ND*ND+j*ND+k]=c2h[i*ND*ND+j*ND+k]+c2ddtt*delt;//Time evolution of the concentration field
				c3h2[i*ND*ND+j*ND+k]=c3h[i*ND*ND+j*ND+k]+c3ddtt*delt;

				//if(c2h[i][j][k]>=1.0){c2h[i][j][k]=1.0-1.0e-06;}//Domain correction of the concentration field
				//if(c2h[i][j][k]<=0.0){c2h[i][j][k]=1.0e-06;}
				//if(c3h[i][j][k]>=1.0){c3h[i][j][k]=1.0-1.0e-06;}
				//if(c3h[i][j][k]<=0.0){c3h[i][j][k]=1.0e-06;}
				if(c2h[i*ND*ND+j*ND+k]>=1.0){c2h[i*ND*ND+j*ND+k]=1.0-1.0e-06;}//Domain correction of the concentration field
				if(c2h[i*ND*ND+j*ND+k]<=0.0){c2h[i*ND*ND+j*ND+k]=1.0e-06;}
				if(c3h[i*ND*ND+j*ND+k]>=1.0){c3h[i*ND*ND+j*ND+k]=1.0-1.0e-06;}
				if(c3h[i*ND*ND+j*ND+k]<=0.0){c3h[i*ND*ND+j*ND+k]=1.0e-06;}
			}
		}
	}

//*** Concentration field balance correction ***********************************************
	//sumc2=0.; for(i=0;i<=ndm;i++){  for(j=0;j<=ndm;j++){ sumc2+=c2h2[i][j][k]; } }
	//dc2a=sumc2/ND/ND/ND-c2a;
	//for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ for(k=0;k<=ndm;k++){ c2h[i][j][k]=c2h2[i][j][k]-dc2a; } } }
	sumc2=0.;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				sumc2+=c2h2[i*ND*ND+j*ND+k];
			}
		}
	}
	dc2a=sumc2/ND/ND/ND-c2a;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				c2h[i*ND*ND+j*ND+k]=c2h2[i*ND*ND+j*ND+k]-dc2a;
			}
		}
	}

	//sumc3=0.; for(i=0;i<=ndm;i++){  for(j=0;j<=ndm;j++){ for(k=0;k<=ndm;k++){ sumc3+=c3h2[i][j][k]; } } }
	//dc3a=sumc3/ND/ND/ND-c3a;
	//for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ for(k=0;k<=ndm;k++){ c3h[i][j][k]=c3h2[i][j][k]-dc3a; } } }
	sumc3=0.;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				sumc3+=c3h2[i*ND*ND+j*ND+k];
			}
		}
	}
	dc3a=sumc3/ND/ND/ND-c3a;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				c3h[i*ND*ND+j*ND+k]=c3h2[i*ND*ND+j*ND+k]-dc3a;
			}
		}
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//sumc23=c2h[i][j][k]+c3h[i][j][k];
				//if(sumc23>=1.0){
				//	c2h[i][j][k]=c2h[i][j][k]/sumc23-1.0e-06;
				//	c3h[i][j][k]=c3h[i][j][k]/sumc23-1.0e-06;
				//}
				sumc23=c2h[i*ND*ND+j*ND+k]+c3h[i*ND*ND+j*ND+k];
				if(sumc23>=1.0){
					c2h[i*ND*ND+j*ND+k]=c2h[i*ND*ND+j*ND+k]/sumc23-1.0e-06;
					c3h[i*ND*ND+j*ND+k]=c3h[i*ND*ND+j*ND+k]/sumc23-1.0e-06;
				}
			}
		}
	}
//*********************************************************************

	//if(keypress()){return 0;}	//Waiting for key
	time1=time1+1.0;
	
	if(time1<time1max){goto start;}//Determining if the maximum count has been reached
	printf("Finished \n");

end:;
	std::exit(0);
	//return 0;
}


//************ Initial concentration field setting subroutine *************
void ini000(double *c2h, double *c3h, int ND)
{
	int i, j, k, id;
 	//srand(time(NULL));//random number seeding
	int ndm=ND-1;

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//c2h[i][j][k]=c2a+0.01*(2.0*DRND(1)-1.0);//Set the concentration field with random numbers up to }1%
				//c3h[i][j][k]=c3a+0.01*(2.0*DRND(1)-1.0);
				c2h[i*ND*ND+j*ND+k]=c2a+0.01*(2.0*DRND(1)-1.0);//Set the concentration field with random numbers up to }1%
				c3h[i*ND*ND+j*ND+k]=c3a+0.01*(2.0*DRND(1)-1.0);
			}
		}
	}

}

//************ data save subroutine *******************************
void datsave(double *c2h, double *c3h, int ND)
{
	FILE *stream;	//Stream pointer setting
	char	fName[256];
	int i, j, k;	//integer
	int ndm=ND-1;

	sprintf(fName,"data_%06d.dat",iout);
	//stream = fopen("test.dat", "a");	//Open the file to write to in append mode
	stream = fopen(fName, "w");
	fprintf(stream, "%d %d %d \n", ndm, ndm, ndm);
	fprintf(stream, "%f\n", time1);		//Saving calculation counts
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//fprintf(stream, "%e  %e  ", c2h[i][j][k], c3h[i][j][k]);//Conservation of local concentration fields
				fprintf(stream, "%e  %e  ", c2h[i*ND*ND+j*ND+k], c3h[i*ND*ND+j*ND+k]);//Conservation of local concentration fields
			}
		}
	}
	fprintf(stream, "\n");	//writing a newline
	fclose(stream);					//close file
}

void datsave_paraview(double *c2h, double *c3h, int ND)
{
	FILE	*fp;
	char	fName[256];
	int 	i, j, k;
	int ndm=ND-1;
	
	printf("sp_result%06d.vtk \n",iout);
	sprintf(fName,"sp_result%06d.vtk",iout);
	fp = fopen(fName, "w");
	fprintf(fp,"# vtk DataFile Version 3.0 \n");
	fprintf(fp,"output.vtk \n");
	fprintf(fp,"ASCII \n");
	fprintf(fp,"DATASET STRUCTURED_POINTS \n");
	fprintf(fp,"DIMENSIONS %5d %5d %5d \n",(ndm+1),(ndm+1),(ndm+1));
	fprintf(fp,"ORIGIN 0.0 0.0 0.0 \n");
	fprintf(fp,"ASPECT_RATIO 1 1 1 \n");
	fprintf(fp,"POINT_DATA %16d \n",((ndm+1)*(ndm+1)*(ndm+1)));
	fprintf(fp,"SCALARS concentration_A float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndm;k++){
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				//fprintf(stream, "%e  %e  ", c2h[i][j][k], c3h[i][j][k]);//Conservation of local concentration fields
				//fprintf(fp,"%10.6f\n", (1.0-c2h[i][j][k]-c3h[i][j][k]));
				fprintf(fp,"%10.6f\n", (1.0-c2h[i*ND*ND+j*ND+k]-c3h[i*ND*ND+j*ND+k]));
			}
		}
	}
	fprintf(fp,"SCALARS concentration_B float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndm;k++){
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				//fprintf(stream, "%e  %e  ", c2h[i][j][k], c3h[i][j][k]);//Conservation of local concentration fields
				//fprintf(fp,"%10.6f\n", c2h[i][j][k]);
				fprintf(fp,"%10.6f\n", c2h[i*ND*ND+j*ND+k]);
			}
		}
	}
	fprintf(fp,"SCALARS concentration_C float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndm;k++){
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				//fprintf(stream, "%e  %e  ", c2h[i][j][k], c3h[i][j][k]);//Conservation of local concentration fields
				//fprintf(fp,"%10.6f\n", c3h[i][j][k]);
				fprintf(fp,"%10.6f\n", c3h[i*ND*ND+j*ND+k]);
			}
		}
	}
	fclose(fp);
}
//************ Reading field data *******************************
void datin(double *c2h, double *c3h, int ND)
{
	FILE *datin0;//Stream pointer setting
	int i, j, k;	//integer
	int ndm=ND-1;
	int rndxm, rndym, rndzm;

	datin0 = fopen("data.dat", "r");//Open the source file

	fscanf(datin0, "%d %d %d", &rndxm, &rndym, &rndzm);
	if (ndm != rndxm){
		printf("data size is mismatch \n");
		printf("Please, change ND= %d in parameters.txt \n", rndxm);
	}
	
	fscanf(datin0, "%lf", &time1);
	
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				fscanf(datin0, "%e  %e  ", &c2h[i*ND*ND+j*ND+k], &c3h[i*ND*ND+j*ND+k]);//Conservation of local concentration fields
			}
		}
	}
	printf("time=  %f  \n", time1);
	fclose(datin0);			//close file
}