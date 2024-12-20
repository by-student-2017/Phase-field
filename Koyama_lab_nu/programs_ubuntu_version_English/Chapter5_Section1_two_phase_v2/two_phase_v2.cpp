#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

#define DRND(x) ((double)(x)/RAND_MAX*rand())	//Random number setting

//#define ND 128				//Number of divisions per side of calculation domain in difference calculation

	//int nd=ND;				//Number of difference divisions on one side of the computational domain (number of difference blocks)
	//int ndm=ND-1;				/define /ND-1
	double PI=3.141592654, rr=8.3145, temp, time1;//Pi, gas constant, absolute temperature, number of calculation counts (proportional to time)
	double c2a;					//Average composition (B component)
	//double c2h[ND][ND];		//topical composition
	//double ph[ND][ND];		//local volume fraction
	int Nstep, iout;

	void ini000(double *c2h, int ND);				//Initial concentration profile setting subroutine
	void datsave(double *ph, double *c2h, int ND);				//data save subroutine
	void datsave_paraview(double *ph, double *c2h, int ND);	//data save subroutine

//******* main program ******************************************
int main(void)
{
	int ND, nd, ndm;
	
	double c1, c2;				//local concentration
	double c2max, c2min;		//Maximum and minimum density
	//double c2h2[ND][ND];		//auxiliary matrix for the local concentration field
	//double c2k_chem, c2k_su, c2k[ND][ND];	//local potential
	double c2k_chem, c2k_su;	//local potential
	double dc2a, sumc2;			//Variables used for concentration field balance calculation
	double dakd2;				//Second derivative of the diffusion potential
	double c2ddtt;				//Time variation of concentration field

	int   i, j;					//integer
	int   ip, im, jp, jm;		//(i+1),(i-1),(j+1),(j-1)
	double al, b1, rtemp, delt;	//Length of one side of computational domain, length of one side of difference block, RT, time step
	double time1max;			//Maximum calculation count (calculation end count)
	double cmob22;				//mobility

	double om_12, om_12e;		//interaction parameter
	double kapa_c2, kapa_c2c;	//concentration gradient energy coefficient

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
	c2a     = data[1];
	delt    = data[2];	// timestep
	temp    = data[3];	// [K]
	al      = data[4];	// [nm]
	cmob22  = data[5];
	om_12e  = data[6];
	kapa_c2c= data[7];
	time1max= int(data[8]);
	Nstep   = int(data[9]);
	printf("---------------------------------\n");
	//
	nd=ND;
	ndm=ND-1;
	//
	//double ch[(nstep+1)*nstep]; // [(nstep+1)*nstep]=[nstep*nstep+nstep], a[i][j]= a[i*n+j] for a[][n]
	//a[z * ySize * xSize + y * xSize + x], a[i][j][k]=a[i*n*m + j*n + k]
	double *c2h  = (double *)malloc(sizeof(double)*( ND*ND ));	//topical composition
	double *ph   = (double *)malloc(sizeof(double)*( ND*ND ));	//local volume fraction
	double *c2h2 = (double *)malloc(sizeof(double)*( ND*ND ));	//auxiliary matrix for the local concentration field
	double *c2k  = (double *)malloc(sizeof(double)*( ND*ND ));	//local potential

	//printf("C2B(0.45) =  ");	scanf(" %lf",&c2a);//Input average composition (B component) from standard input/output	//	c2a=0.45;
	//printf("delt(0.005)=  ");	scanf(" %lf",&delt);//ticking	//	delt=0.005;

	//temp=900.0;				//temperature (K)
	rtemp=rr*temp;				//RT
	//al=100.0*1.0E-09;			//Side length of calculation area (nm)
	al=al*1.0E-09;				//Side length of calculation area (nm)
	b1=al/(double)ND;			//Length of one side of difference block

	//cmob22=1.0;				//mobility

	//om_12=25000./rtemp; 		//Interaction parameters (in J/mol, non-dimensionalized at RT)
	om_12=om_12e/rtemp; 		//Interaction parameters (in J/mol, non-dimensionalized at RT)

	//kapa_c2=5.0e-15/b1/b1/rtemp;//gradient energy coefficient (in Jm^2/mol, non-dimensionalized by RT and b1^2)
	kapa_c2=kapa_c2c/b1/b1/rtemp;//gradient energy coefficient (in Jm^2/mol, non-dimensionalized by RT and b1^2)

	time1=0.0;
	//time1max=1.0e05+1.0;//Initial and maximum calculation counts

//*** Initial concentration field setting and drawing window display *****************************************

	ini000(c2h, ND);					//Setting the initial concentration field

//**** Simulation start ******************************
//Nstep = 1000;
iout = -1;
start: ;

	if((((int)(time1) % Nstep)==0)) {datsave(ph, c2h, ND);} //Save the concentration field every fixed repetition count
	if((((int)(time1) % Nstep)==0)) {datsave_paraview(ph, c2h, ND);} //Save the concentration field every fixed repetition count

//***** Computation of the potential field (see Section 5.2) ***********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ip=i+1; im=i-1; jp=j+1; jm=j-1;
			if(i==ndm){ip=0;}  if(i==0){im=ndm;}//periodic boundary conditions
			if(j==ndm){jp=0;}  if(j==0){jm=ndm;}

			//c2=c2h[i][j];  c1=1.0-c2;//local concentration field
			c2=c2h[i*ND+j];  c1=1.0-c2;//local concentration field

			c2k_chem=om_12*(c1-c2)+(log(c2)-log(c1));//chemical diffusion potential
			//c2k_su=-2.*kapa_c2*(c2h[ip][j]+c2h[im][j]+c2h[i][jp]+c2h[i][jm]-4.0*c2);//gradient potential
			c2k_su=-2.*kapa_c2*(c2h[ip*ND+j]+c2h[im*ND+j]+c2h[i*ND+jp]+c2h[i*ND+jm]-4.0*c2);//gradient potential

			//c2k[i][j]=c2k_chem+c2k_su;//diffusion potential
			c2k[i*ND+j]=c2k_chem+c2k_su;//diffusion potential
		}
	}

//***** Evolution equation calculation (see section 5.2) **********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ip=i+1; im=i-1; jp=j+1; jm=j-1;
			if(i==ndm){ip=0;} if(i==0){im=ndm;}//periodic boundary conditions
			if(j==ndm){jp=0;} if(j==0){jm=ndm;}

			//dakd2=c2k[ip][j]+c2k[im][j]+c2k[i][jp]+c2k[i][jm]-4.0*c2k[i][j];//Second derivative of the diffusion potential
			dakd2=c2k[ip*ND+j]+c2k[im*ND+j]+c2k[i*ND+jp]+c2k[i*ND+jm]-4.0*c2k[i*ND+j];//Second derivative of the diffusion potential

			c2ddtt=cmob22*dakd2;//diffusion equation

			//c2h2[i][j]=c2h[i][j]+c2ddtt*delt;//Time evolution of the concentration field
			c2h2[i*ND+j]=c2h[i*ND+j]+c2ddtt*delt;//Time evolution of the concentration field

			//if(c2h[i][j]>=1.0){c2h[i][j]=1.0-1.0e-06;}//Concentration field correction
			//if(c2h[i][j]<=0.0){c2h[i][j]=1.0e-06;}
			if(c2h[i*ND+j]>=1.0){c2h[i*ND+j]=1.0-1.0e-06;}//Concentration field correction
			if(c2h[i*ND+j]<=0.0){c2h[i*ND+j]=1.0e-06;}
		}
	}

//*** Concentration field balance correction ***********************************************
	//sumc2=0.; for(i=0;i<=ndm;i++){  for(j=0;j<=ndm;j++){ sumc2+=c2h2[i][j]; } }
	sumc2=0.; for(i=0;i<=ndm;i++){  for(j=0;j<=ndm;j++){ sumc2+=c2h2[i*ND+j]; } }
	dc2a=sumc2/ND/ND-c2a;
	//for(i=0;i<=ndm;i++){  for(j=0;j<=ndm;j++){ c2h[i][j]=c2h2[i][j]-dc2a; } }
	for(i=0;i<=ndm;i++){  for(j=0;j<=ndm;j++){ c2h[i*ND+j]=c2h2[i*ND+j]-dc2a; } }

//*** Setting the phase field (see section 5.3.1) ***************************************
	c2max=0.0; c2min=1.0;
	for(i=0;i<=ndm;i++){  
		for(j=0;j<=ndm;j++){ 
			//if(c2max<=c2h[i][j]){c2max=c2h[i][j];}
			//if(c2min>=c2h[i][j]){c2min=c2h[i][j];}
			if(c2max<=c2h[i*ND+j]){c2max=c2h[i*ND+j];}
			if(c2min>=c2h[i*ND+j]){c2min=c2h[i*ND+j];}
		} 
	}

	for(i=0;i<=ndm;i++){  
		for(j=0;j<=ndm;j++){ 
			//ph[i][j]=(c2h[i][j]-c2min)/(c2max-c2min);
			ph[i*ND+j]=(c2h[i*ND+j]-c2min)/(c2max-c2min);
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
void ini000(double *c2h, int ND)
{
	int i, j, id;
 	//srand(time(NULL));//random number seeding
	int ndm=ND-1;

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//c2h[i][j]=c2a+0.01*(2.0*DRND(1)-1.0);//Set the concentration field with random numbers up to �}1%
			c2h[i*ND+j]=c2a+0.01*(2.0*DRND(1)-1.0);
		}
	}

}

//************ data save subroutine *******************************
void datsave(double *ph, double *c2h, int ND)
{
	FILE *stream;	//Stream pointer setting
	int i, j;			//integer
	int ndm=ND-1;

	stream = fopen("ph.dat", "a");	//Open the file to write to in append mode
	fprintf(stream, "%f\n", time1);		//Save repeat count
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//fprintf(stream, "%e  ", ph[i][j]);//Save phase field
			fprintf(stream, "%e  ", ph[i*ND+j]);
			//fprintf(stream, "%e  ", c2h[i][j]);//Conservation of local concentration fields
		}
	}
	fprintf(stream, "\n");	//writing a newline
	fclose(stream);					//close file
}

void datsave_paraview(double *ph, double *c2h, int ND)
{
	FILE	*fp;
	char	fName[256];
	int 	i, j, k;
	int ndm=ND-1;
	
	iout = iout + 1;
	printf("pf_result%06d.vtk \n",iout);
	sprintf(fName,"pf_result%06d.vtk",iout);
	fp = fopen(fName, "w");
	fprintf(fp,"# vtk DataFile Version 3.0 \n");
	fprintf(fp,"output.vtk \n");
	fprintf(fp,"ASCII \n");
	fprintf(fp,"DATASET STRUCTURED_POINTS \n");
	fprintf(fp,"DIMENSIONS %5d %5d %5d \n",(ndm+1),(ndm+1),1);
	fprintf(fp,"ORIGIN 0.0 0.0 0.0 \n");
	fprintf(fp,"ASPECT_RATIO 1 1 1 \n");
	fprintf(fp,"POINT_DATA %11d \n",((ndm+1)*(ndm+1)*1));
	fprintf(fp,"SCALARS Phase_field float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			//fprintf(fp, "%e  ", ph[i][j]);//Save phase field
			//fprintf(fp,"%10.6f\n", ph[i][j]);
			fprintf(fp,"%10.6f\n", ph[i*ND+j]);
		}
	}
	fprintf(fp,"SCALARS concentration float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			//fprintf(fp, "%e  ", c2h[i][j]);//Conservation of local concentration fields
			//fprintf(fp,"%10.6f\n", c2h[i][j]);
			fprintf(fp,"%10.6f\n", c2h[i*ND+j]);
		}
	}
	fclose(fp);
}