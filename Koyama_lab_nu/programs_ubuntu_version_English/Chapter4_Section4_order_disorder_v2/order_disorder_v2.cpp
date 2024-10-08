#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

#define DRND(x) ((double)(x)/RAND_MAX*rand())//Random number function setting

//#define ND 400					//Number of divisions per side of calculation domain in difference calculation

	//int nd=ND, ndm=ND-1; 			//Define the number of difference divisions (number of difference blocks) on one side of the calculation area, ND-1
	//int nd2=ND/2;				 	//Defined use ND/2
	double PI=3.14159;				//Pi
	double rr=8.3145;				//gas constant
	double time1;					//time
	int iout;

	//double s1h[ND][ND];				//long range regular field

	void ini000(double *s1h, int ND);	//Initial field setup subroutine
	void datsave(double *s1h, int ND);	//data save subroutine
	void datsave_paraview(double *s1h, int ND);//data save subroutine

//******* main program ******************************************
int main(void)
{
	int ND, nd, ndm, nd2;
	
	double s1;						//long range regular field
	//double s1h2[ND][ND];			//Auxiliary array of fields
	double s1k_chem, s1k_su;		//potential
	double sum11;					//spatial integral of s1
	double s1ddtt;					//Time variation of s1 (left side of evolution equation)

	int   i, j, k, l, ii, jj, kk, iii, jjj;	//integer
	int   p, q, m, n;						//integer
	int   ip, im, jp, jm, Nstep;			//integer
	double al, temp, delt;					//Computational domain, time, time step
	double time1max;						//Maximum time (used to stop calculation)
	double b1, vm0, atom_n;					//normalized length, molar volume, number of atoms in unit cell
	double smob;							//Relaxation coefficient of crystal transformation

	double AA0, AA1, AA2, AA3;				//Coefficient in Gisb Energy
	double AA0e;
	double a1_c, b1_c, c1_c;				//lattice constant
	double kappa_s1;						//gradient energy factor
	double kappa_s1c;
	//double ds_fac;						//Fluctuation coefficient of crystal transformation

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
	delt    = data[1];
	temp    = data[2];
	al      = data[3];	// [nm]
	smob    = data[4];
	AA0e    = data[5];
	kappa_s1c=data[6];
	vm0     = data[7];
	time1max= int(data[8]);
	Nstep   = int(data[9]);
	printf("---------------------------------\n");
	
	//
	nd=ND, ndm=ND-1; 			//Define the number of difference divisions (number of difference blocks) on one side of the calculation area, ND-1
	nd2=ND/2;				 	//Defined use ND/2
	//
	double *s1h  = (double *)malloc(sizeof(double)*( ND*ND ));	//long range regular field
	double *s1h2 = (double *)malloc(sizeof(double)*( ND*ND ));	//Auxiliary array of fields
	//
	//printf("DELT(1.0)=  ");	scanf(" %lf",&delt);//Time step input	//delt=1.0;

	//temp=500.0;					//temperature (K)
	//al=1000.0*1.0E-09;			//Computational domain [m]
	al=al*1.0E-09;				//Computational domain [m]
	b1=al/nd;					//Diff block length

	time1=0.0;					//Initial calculation count setting
	//time1max=1.0+1.0e+07;		//Setting the maximum calculation count

	//smob=1.0;					//Mobility (relaxation coefficient of transformation)

	//AA0=200.0/rr/temp;			//The chemical driving force for rule-disorder transformations
	AA0=AA0e/rr/temp;			//The chemical driving force for rule-disorder transformations

	//kappa_s1=5.0e-15/rr/temp/b1/b1;//gradient energy factor
	kappa_s1=kappa_s1c/rr/temp/b1/b1;//gradient energy factor

	//a1_c=b1_c=c1_c=3.563E-10;	//lattice constant
	//atom_n=4.0;  vm0=6.02E23*a1_c*b1_c*c1_c/atom_n;	//Calculation of molar volume (assuming fcc)

//*** Initial field setting ***************

	ini000(s1h, ND);

//**** Simulation start ******************************
//Nstep = 10;
iout = -1;
start: ;

	//if(time1<=200.){Nstep=10;} else{Nstep=100;}		//Changing the time interval for saving data
	if((((int)(time1) % Nstep)==0)) {datsave(s1h, ND);} 	//Save tissue data every fixed repeat count
	if((((int)(time1) % Nstep)==0)) {datsave_paraview(s1h, ND);} 	//Save tissue data every fixed repeat count

//******  Potential calculation ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ip=i+1; im=i-1; jp=j+1; jm=j-1;
			if(i==ndm){ip=0;}  if(i==0){im=ndm;}
			if(j==ndm){jp=0;}  if(j==0){jm=ndm;}

			//s1=s1h[i][j];
			//s1k_su=-kappa_s1*(s1h[ip][j]+s1h[im][j]+s1h[i][jp]+s1h[i][jm]-4.0*s1);
			s1=s1h[i*ND+j];
			s1k_su=-kappa_s1*(s1h[ip*ND+j]+s1h[im*ND+j]+s1h[i*ND+jp]+s1h[i*ND+jm]-4.0*s1);
															//gradient potential [equation (4.16)]
			s1k_chem=-4.0*AA0*s1*(1.0-s1*s1);				//Calculation of chemical potential [equation (4.14)]

			s1ddtt=-smob*(s1k_chem+s1k_su);					//Calculation of the time evolution of the field [equation (4.17)]
			//s1h2[i][j]=s1h[i][j]+s1ddtt*delt;				//Explicit method
			s1h2[i*ND+j]=s1h[i*ND+j]+s1ddtt*delt;			//Explicit method

			//if(s1h2[i][j]>=1.0){s1h2[i][j]=1.0;}  
			//if(s1h2[i][j]<=-1.0){s1h2[i][j]=-1.0;}			//Correction for the domain of s (-1<=s<=1)
			if(s1h2[i*ND+j]>= 1.0){s1h2[i*ND+j]= 1.0;}  
			if(s1h2[i*ND+j]<=-1.0){s1h2[i*ND+j]=-1.0;}			//Correction for the domain of s (-1<=s<=1)
		}
	}

	//for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){  s1h[i][j]=s1h2[i][j]; } }
	for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){  s1h[i*ND+j]=s1h2[i*ND+j]; } }
	//Copy data from auxiliary array to main array

	//if(keypress()){return 0;}				//Waiting for key

	time1=time1+1.0;								//Add calculation count
	if(time1<time1max){goto start;}	//Determining if the maximum count has been reached
	printf("Finished \n");

end:;
	std::exit(0);
	//return 0;
}

//************ Initial wave setting subroutine *******************************
void ini000(double *s1h, int ND)
{
	int i, j;
	srand(time(NULL)); // random number initialization
	int ndm=ND-1;

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//s1h[i][j]=0.01*(2.0*DRND(1.0)-1.0);//Set the field with random numbers up to �}1%
			s1h[i*ND+j]=0.01*(2.0*DRND(1.0)-1.0);//Set the field with random numbers up to �}1%
		}
	}
}

//************ Save data subroutine *******************************
void datsave(double *s1h, int ND)
{
	FILE		*stream;	//Stream pointer setting
	int 		i, j;			//integer
	int ndm=ND-1;

	stream = fopen("test.dat", "a");	//Open the write destination file for appending

	fprintf(stream, "%e\n", time1);		//Save repeat count
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//fprintf(stream, "%e  ", s1h[i][j]);//Field data storage
			fprintf(stream, "%e  ", s1h[i*ND+j]);//Field data storage
		}
	}
	fprintf(stream, "\n");	//writing a newline
	fclose(stream);					//close file
}

void datsave_paraview(double *s1h, int ND)
{
	FILE	*fp;
	char	fName[256];
	int 	i, j;
	int ndm=ND-1;
	
	iout = iout + 1;
	printf("pf_result%06d.vtk \n",iout);
	sprintf(fName,"pf_result%06d.vtk",iout);
	fp = fopen(fName, "w");
	fprintf(fp,"# vtk DataFile Version 3.0 \n");
	fprintf(fp,"output.vtk \n");
	fprintf(fp,"ASCII \n");
	fprintf(fp,"DATASET STRUCTURED_POINTS \n");
	fprintf(fp,"DIMENSIONS %5d %5d %5d \n",(ndm+1),(ndm+1),(1));
	fprintf(fp,"ORIGIN 0.0 0.0 0.0 \n");
	fprintf(fp,"ASPECT_RATIO 1 1 1 \n");
	fprintf(fp,"POINT_DATA %16d \n",((ndm+1)*(ndm+1)*1));
	fprintf(fp,"SCALARS Phase_field float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			//fprintf(stream, "%e  ", s1h[i][j]);//Field data storage
			//fprintf(fp,"%10.6f\n", s1h[i][j]);
			fprintf(fp,"%10.6f\n", s1h[i*ND+j]);
		}
	}
	fclose(fp);
}
