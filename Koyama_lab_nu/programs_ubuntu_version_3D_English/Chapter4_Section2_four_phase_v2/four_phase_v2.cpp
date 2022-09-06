#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

#define DRND(x) ((double)(x)/RAND_MAX*rand())  //Random number setting

	double PI=3.141592654;		//Pi
	double rr=8.3145;			//gas constant
	double temp;				//absolute temperature
	double time1;				//Number of calculation counts (proportional to time)
	double c2a, c3a, c4a;		//Average composition (1: A component, 2: B component, 3: C component)
	int Nstep, iout;

	void ini000(double *c2h, double *c3h, double *c4h, int ND);	//Initial concentration profile setting subroutine
	void datin(double *c2h, double *c3h, double *c4h, int ND);	//Subroutine for initial field reading
	void datsave(double *c2h, double *c3h, double *c4h, int ND);	//data save subroutine
	void datsave_paraview(double *c2h, double *c3h, double *c4h, int ND);//data save subroutine

//******* main program ******************************************
int main(void)
{
	int ND, nd, ndm;
	
	double c1, c2, c3, c4;							//local concentration
	double c2k_chem, c2k_su;					//local potential
	double c3k_chem, c3k_su;					//local potential
	double c4k_chem, c4k_su;					//local potential
	double dc2a, sumc2;	//Variables used for concentration field balance calculation
	double dc3a, sumc3;
	double dc4a, sumc4;
	double sumct;
	double dakd2, dakd3, dakd4;					//Second derivative of the diffusion potential
	double c2ddtt, c3ddtt, c4ddtt;				//Time variation of concentration field

	int   i, j, k;								//integer
	int   ip, im, jp, jm, kp, km;				//(i+1),(i-1),(j+1),(j-1),(k+1),(k-1)
	double al, b1, rtemp, delt;					//Length of one side of computational domain, length of one side of difference block, RT, time step
	double time1max;							//Maximum calculation count (calculation end count)
	double cmob22, cmob33, cmob44;				//mobility
	double cmob23, cmob32;
	double cmob24, cmob42;
	double cmob34, cmob43;

	double om_12, om_13, om_14;					//interaction parameter
	double om_23, om_24, om_34;
	double om_12e, om_13e, om_14e;				//interaction parameter
	double om_23e, om_24e, om_34e;
	double kapa_c1, kapa_c2, kapa_c3, kapa_c4;		//concentration gradient energy coefficient
	double kapa_c1c, kapa_c2c, kapa_c3c, kapa_c4c;	//concentration gradient energy coefficient
	
	double div2_c2h, div2_c3h, div2_c4h;		//div(div(ch))=d^2(ch)/dx^2 + d^2(ch)/dy^2 + d^2(ch)/dz^2
	
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
	c4a     = data[3];
	delt    = data[4];
	temp    = data[5];
	al      = data[6];	// [nm]
	cmob22  = data[7];
	cmob33  = data[8];
	cmob44  = data[9];
	cmob23  = data[10];
	cmob24  = data[11];
	cmob34  = data[12];
	om_12e  = data[13];	//L_AB
	om_13e  = data[14];	//L_AC
	om_14e  = data[15];	//L_AD
	om_23e  = data[16];	//L_BC
	om_24e  = data[17];	//L_BD
	om_34e  = data[18];	//L_CD
	kapa_c1c= data[19];
	kapa_c2c= data[20];
	kapa_c3c= data[21];
	kapa_c4c= data[22];
	time1max= int(data[23]);
	Nstep   = int(data[24]);
	readff  = int(data[25]);
	printf("---------------------------------\n");
	//
	nd=ND;					//Number of difference divisions on one side of the computational domain (number of difference blocks)
	ndm=ND-1;				//Define ND-1
	//
	//double ch[(nstep+1)*nstep]; // [(nstep+1)*nstep]=[nstep*nstep+nstep], a[i][j]= a[i*n+j] for a[][n]
	//a[z * ySize * xSize + y * xSize + x], a[i][j][k]=a[i*n*m + j*n + k]
	double *c2h  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//topical composition
	double *c3h  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//topical composition
	double *c4h  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//topical composition
	double *c2h2 = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//auxiliary matrix for the local concentration field
	double *c3h2 = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//auxiliary matrix for the local concentration field
	double *c4h2 = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//auxiliary matrix for the local concentration field
	double *c2k  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//local potential
	double *c3k  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//local potential
	double *c4k  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//local potential
	//
	rtemp=rr*temp;				//RT
	
	al=al*1.0E-09;				//Side length of calculation area (m)
	b1=al/(double)ND;			//Length of one side of difference block

	cmob32=cmob23;				//mobility
	cmob42=cmob24;
	cmob34=cmob43;

	om_12=om_12e/rtemp; 		//Interaction parameters (in J/mol, non-dimensionalized at RT)
	om_13=om_13e/rtemp;
	om_14=om_14e/rtemp;
	om_23=om_23e/rtemp;
	om_24=om_24e/rtemp;
	om_34=om_34e/rtemp;

	kapa_c1=kapa_c1c/b1/b1/rtemp;//gradient energy coefficient (in Jm^2/mol, non-dimensionalized by RT and b1^2)
	kapa_c2=kapa_c2c/b1/b1/rtemp;//gradient energy coefficient (in Jm^2/mol, non-dimensionalized by RT and b1^2)
	kapa_c3=kapa_c3c/b1/b1/rtemp;//gradient energy coefficient (in Jm^2/mol, non-dimensionalized by RT and b1^2)
	kapa_c4=kapa_c4c/b1/b1/rtemp;//gradient energy coefficient (in Jm^2/mol, non-dimensionalized by RT and b1^2)

	time1=0.0;					//Initial value of calculation count

//*** Initial concentration field setting and drawing window display *****************************************
	if(readff==0){
		printf("make initial concentration (random) \n");
		ini000(c2h, c3h, c4h, ND);//Setting the initial concentration field
	} else {
		printf("read data.dat file \n");
		datin(c2h, c3h, c4h, ND);
	}

//**** Simulation start ******************************
//Nstep = 2000;
iout = -1;
start: ;

	//printf("Time: %f \n", time1);
	if((((int)(time1) % Nstep)==0)) {iout = iout + 1;}
	if((((int)(time1) % Nstep)==0)) {datsave(c2h, c3h, c4h, ND);} //Save the concentration field every fixed repetition count
	if((((int)(time1) % Nstep)==0)) {datsave_paraview(c2h, c3h, c4h, ND);} //Save the concentration field every fixed repetition count

//***** Potential field calculation ***********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				
				ip=i+1; im=i-1; jp=j+1; jm=j-1; kp=k+1; km=k-1;
				if(i==ndm){ip=0;} if(i==0){im=ndm;}	//periodic boundary conditions
				if(j==ndm){jp=0;} if(j==0){jm=ndm;}	//periodic boundary conditions
				if(k==ndm){kp=0;} if(k==0){km=ndm;}	//periodic boundary conditions

				//c2=c2h[i][j]; c3=c3h[i][j]; c1=1.0-c2-c3;//local concentration field
				c2=c2h[i*ND*ND+j*ND+k]; c3=c3h[i*ND*ND+j*ND+k]; c4=c4h[i*ND*ND+j*ND+k]; 
				c1=1.0-c2-c3-c4;//local concentration field

				//Gsys = Gc + Ggrad
				//Gc = integral { (L_AB * c1 * c2) + (L_AC * c1 * c3) + (L_BC * c2 * c3) + ...
				//  + RT*( c1*log(c1) + c2*log(c2) + c3*log(c3) + ...) }  dr
				c2k_chem=om_12*c1-om_12*c2-om_13*c3-om_14*c4+om_23*c3+om_24*c4+(log(c2)-log(c1));//chemical diffusion potential

				div2_c2h = (c2h[ip*ND*ND+j*ND+k]+c2h[im*ND*ND+j*ND+k]
						   +c2h[i*ND*ND+jp*ND+k]+c2h[i*ND*ND+jm*ND+k]
      					   +c2h[i*ND*ND+j*ND+kp]+c2h[i*ND*ND+j*ND+km]
						   -6.0*c2);
				
				div2_c3h = (c3h[ip*ND*ND+j*ND+k]+c3h[im*ND*ND+j*ND+k]
						   +c3h[i*ND*ND+jp*ND+k]+c3h[i*ND*ND+jm*ND+k]
      					   +c3h[i*ND*ND+j*ND+kp]+c3h[i*ND*ND+j*ND+km]
						   -6.0*c3);
								
				div2_c4h = (c4h[ip*ND*ND+j*ND+k]+c4h[im*ND*ND+j*ND+k]
						   +c4h[i*ND*ND+jp*ND+k]+c4h[i*ND*ND+jm*ND+k]
      					   +c4h[i*ND*ND+j*ND+kp]+c4h[i*ND*ND+j*ND+km]
						   -6.0*c4);
				
				//0.5*(c2+c3+c4)*(c2+c3+c4)+0.5*c2*c2+c3*c3+c4*c4
				//						 =0.5*c2*c2+c2*c3+c2*c4
				//						 +0.5*c3*c2+c3*c3+c3+c4
				//						 +0.5*c4*c2+c4*c3+c4+c4
				//						 +0.5*c2*c2+c3*c3+c4*c4
				//						 =c2*c2+c3*c3+c4*c4
				//						 +c2*c3+c2*c4+c3*c4
				//Ggrad = integral {0.5*kappa1*div(div c1) + 0.5*kappa2*div(div c2) + 0.5*kappa3*div(div c3) + ... } dr
				
				//gradient potential
				c2k_su=-2.0*(0.5*kapa_c1+0.5*kapa_c2)*div2_c2h
						   -(0.5*kapa_c1+0.5*kapa_c3)*div2_c3h
						   -(0.5*kapa_c1+0.5*kapa_c4)*div2_c4h;

				c3k_chem=om_13*c1-om_12*c2-om_13*c3-om_14*c4+om_23*c2+om_34*c4+(log(c3)-log(c1));///chemical diffusion potential
				c3k_su=-2.0*(0.5*kapa_c1+0.5*kapa_c3)*div2_c3h
						   -(0.5*kapa_c1+0.5*kapa_c2)*div2_c2h
						   -(0.5*kapa_c1+0.5*kapa_c4)*div2_c4h;
				
				c4k_chem=om_14*c1-om_12*c2-om_13*c3-om_14*c4+om_24*c2+om_34*c3+(log(c4)-log(c1));///chemical diffusion potential
				c4k_su=-2.0*(0.5*kapa_c1+0.5*kapa_c4)*div2_c4h
						   -(0.5*kapa_c1+0.5*kapa_c2)*div2_c2h
						   -(0.5*kapa_c1+0.5*kapa_c3)*div2_c3h;
				
				c2k[i*ND*ND+j*ND+k]=c2k_chem+c2k_su;//Diffusion potential (equation (4.1))
				c3k[i*ND*ND+j*ND+k]=c3k_chem+c3k_su;
				c4k[i*ND*ND+j*ND+k]=c4k_chem+c4k_su;
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

				dakd2=c2k[ip*ND*ND+j*ND+k]+c2k[im*ND*ND+j*ND+k]
					 +c2k[i*ND*ND+jp*ND+k]+c2k[i*ND*ND+jm*ND+k]
      				 +c2k[i*ND*ND+j*ND+kp]+c2k[i*ND*ND+j*ND+km]
					 -6.0*c2k[i*ND*ND+j*ND+k];	//Second derivative of the diffusion potential
				dakd3=c3k[ip*ND*ND+j*ND+k]+c3k[im*ND*ND+j*ND+k]
					 +c3k[i*ND*ND+jp*ND+k]+c3k[i*ND*ND+jm*ND+k]
      				 +c3k[i*ND*ND+j*ND+kp]+c3k[i*ND*ND+j*ND+km]
					 -6.0*c3k[i*ND*ND+j*ND+k];	//Second derivative of the diffusion potential
				dakd4=c4k[ip*ND*ND+j*ND+k]+c4k[im*ND*ND+j*ND+k]
					 +c4k[i*ND*ND+jp*ND+k]+c4k[i*ND*ND+jm*ND+k]
      				 +c4k[i*ND*ND+j*ND+kp]+c4k[i*ND*ND+j*ND+km]
					 -6.0*c4k[i*ND*ND+j*ND+k];	//Second derivative of the diffusion potential

				c2ddtt=cmob22*dakd2+cmob23*dakd3+cmob24*dakd4;//Diffusion equation (equation (4.2))
				c3ddtt=cmob32*dakd2+cmob33*dakd3+cmob34*dakd4;
				c4ddtt=cmob42*dakd2+cmob43*dakd3+cmob44*dakd4;

				c2h2[i*ND*ND+j*ND+k]=c2h[i*ND*ND+j*ND+k]+c2ddtt*delt;//Time evolution of the concentration field
				c3h2[i*ND*ND+j*ND+k]=c3h[i*ND*ND+j*ND+k]+c3ddtt*delt;
				c4h2[i*ND*ND+j*ND+k]=c4h[i*ND*ND+j*ND+k]+c4ddtt*delt;

				if(c2h[i*ND*ND+j*ND+k]>=1.0){c2h[i*ND*ND+j*ND+k]=1.0-1.0e-06;}//Domain correction of the concentration field
				if(c2h[i*ND*ND+j*ND+k]<=0.0){c2h[i*ND*ND+j*ND+k]=1.0e-06;}
				if(c3h[i*ND*ND+j*ND+k]>=1.0){c3h[i*ND*ND+j*ND+k]=1.0-1.0e-06;}
				if(c3h[i*ND*ND+j*ND+k]<=0.0){c3h[i*ND*ND+j*ND+k]=1.0e-06;}
				if(c4h[i*ND*ND+j*ND+k]>=1.0){c4h[i*ND*ND+j*ND+k]=1.0-1.0e-06;}
				if(c4h[i*ND*ND+j*ND+k]<=0.0){c4h[i*ND*ND+j*ND+k]=1.0e-06;}
			}
		}
	}

//*** Concentration field balance correction ***********************************************
	sumc2=0.0; sumc3=0.0; sumc4=0.0;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				sumc2+=c2h2[i*ND*ND+j*ND+k];
				sumc3+=c3h2[i*ND*ND+j*ND+k];
				sumc4+=c4h2[i*ND*ND+j*ND+k];
			}
		}
	}
	dc2a=sumc2/ND/ND/ND-c2a;
	dc3a=sumc3/ND/ND/ND-c3a;
	dc4a=sumc4/ND/ND/ND-c4a;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				c2h[i*ND*ND+j*ND+k]=c2h2[i*ND*ND+j*ND+k]-dc2a;
				c3h[i*ND*ND+j*ND+k]=c3h2[i*ND*ND+j*ND+k]-dc3a;
				c4h[i*ND*ND+j*ND+k]=c4h2[i*ND*ND+j*ND+k]-dc4a;
			}
		}
	}

	sumct=0.0;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				sumct=c2h[i*ND*ND+j*ND+k]+c3h[i*ND*ND+j*ND+k]+c4h[i*ND*ND+j*ND+k];
				if(sumct>=1.0){
					c2h[i*ND*ND+j*ND+k]=c2h[i*ND*ND+j*ND+k]/sumct-1.0e-06;
					c3h[i*ND*ND+j*ND+k]=c3h[i*ND*ND+j*ND+k]/sumct-1.0e-06;
					c4h[i*ND*ND+j*ND+k]=c4h[i*ND*ND+j*ND+k]/sumct-1.0e-06;
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
void ini000(double *c2h, double *c3h, double *c4h, int ND)
{
	int i, j, k, id;
 	//srand(time(NULL));//random number seeding
	int ndm=ND-1;

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				c2h[i*ND*ND+j*ND+k]=c2a+0.01*(2.0*DRND(1)-1.0);//Set the concentration field with random numbers up to }1%
				c3h[i*ND*ND+j*ND+k]=c3a+0.01*(2.0*DRND(1)-1.0);
				c4h[i*ND*ND+j*ND+k]=c4a+0.01*(2.0*DRND(1)-1.0);
			}
		}
	}

}

//************ data save subroutine *******************************
void datsave(double *c2h, double *c3h, double *c4h, int ND)
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
				fprintf(stream, "%e  %e  %e  ", c2h[i*ND*ND+j*ND+k], c3h[i*ND*ND+j*ND+k], c4h[i*ND*ND+j*ND+k]);//Conservation of local concentration fields
			}
		}
	}
	fprintf(stream, "\n");	//writing a newline
	fclose(stream);					//close file
}

void datsave_paraview(double *c2h, double *c3h, double *c4h, int ND)
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
				fprintf(fp,"%10.6f\n", (1.0-c2h[i*ND*ND+j*ND+k]-c3h[i*ND*ND+j*ND+k]-c4h[i*ND*ND+j*ND+k]));
			}
		}
	}
	fprintf(fp,"SCALARS concentration_B float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndm;k++){
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				fprintf(fp,"%10.6f\n", c2h[i*ND*ND+j*ND+k]);
			}
		}
	}
	fprintf(fp,"SCALARS concentration_C float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndm;k++){
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				fprintf(fp,"%10.6f\n", c3h[i*ND*ND+j*ND+k]);
			}
		}
	}
	fprintf(fp,"SCALARS concentration_D float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndm;k++){
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				fprintf(fp,"%10.6f\n", c4h[i*ND*ND+j*ND+k]);
			}
		}
	}
	fprintf(fp,"SCALARS concentration_Total float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndm;k++){
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				fprintf(fp,"%10.6f\n", (1.0-c2h[i*ND*ND+j*ND+k]-c3h[i*ND*ND+j*ND+k]-c4h[i*ND*ND+j*ND+k])*1.0*4.0
									  +c2h[i*ND*ND+j*ND+k]*2.0*4.0
									  +c3h[i*ND*ND+j*ND+k]*3.0*4.0
									  +c4h[i*ND*ND+j*ND+k]*4.0*4.0);
			}
		}
	}
	fclose(fp);
}

//************ Reading field data *******************************
void datin(double *c2h, double *c3h, double *c4h, int ND)
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
				fscanf(datin0, "%le  %le  %le  ", &c2h[i*ND*ND+j*ND+k], &c3h[i*ND*ND+j*ND+k], &c4h[i*ND*ND+j*ND+k]);//Conservation of local concentration fields
			}
		}
	}
	printf("time=  %f  \n", time1);
	fclose(datin0);					//close file
}