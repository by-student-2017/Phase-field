#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

//#include "wingxa.h"

#define DRND(x) ((double)(x)/RAND_MAX*rand())

//#define ND 64
//#define INXY 400					//Pixel size of one side of drawing window

	//int nd=ND, ndm=ND-1, nd2=ND/2; 	//Organization split, organization split-1
	double RR=8.3145;
	double PI=3.141592654;
	double c0, temp, time1;
	//double ch[ND][ND][ND];			//Concentration data array in tissue
	int iout;
	int flg;

	void ini_field(double *ch, int ND);	//Initial concentration wave setting subroutine
	void datin(double *ch, int ND);	//Subroutine for initial field reading
	//void graph_c();					//graph display subroutine
	void datsave(double *ch, int ND);	//Concentration data save subroutine
	void datsave_paraview(double *ch, int ND);	//Concentration data save subroutine
	//void datin();						//Initial concentration wave reading subroutine

//******* Main program ******************************************************
int main(void)
{
	int ND;
	int nd, ndm, nd2; 				//Organization split, organization split-1
	
	//double ck[ND][ND][ND];			//diffusion potential
	//double ch[ND][ND][ND];			//Concentration data array in tissue
	//double ch2[ND][ND][ND];
	double mu_chem, mu_str, mu_surf;
	int    i, j, k, l;
	int    Nstep;
	double time1max; 				//maximum time
	double amob_c;					//atomic mobility constant
	double c;						//concentration
	double cddtt;					//concentration increment
	double kapa_c, kapa_c0;			//concentration gradient energy constant
	double al;						//Computational domain
	double b1;						//Differential block size
	double L0, L00;
	double delt;					//ticking
	double Mx, My, Mz;				//mobility

	double cip, cim, cjp, cjm, ckp, ckm;
	double ck1dev, ck2dev;
	double ck0, ckip, ckim, ckjp, ckjm, ckkp, ckkm;
	int    ip, im, jp, jm, kp, km;
	double sumc, dc0;
	double c_flu, c_flu0;			//Magnitude of concentration field fluctuations
	
	int    readff;

//********************************************************************************
	printf("---------------------------------\n");
	printf("read parameters form parameters.txt\n");
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
	c0      = data[2];
	Mx      = data[3];
	My      = data[4];
	Mz      = data[5];
	al      = data[6];	// [nm]
	amob_c  = data[7];
	time1max= int(data[8]);
	temp    = data[9];	// [K]
	L00     = data[10];	//L0=L00/RR/temp;
	kapa_c0 = data[11];	//kapa_c=kapa_c0/b1/b1/RR/temp;
	c_flu   = data[12];
	flg     = int(data[13]);
	Nstep   = int(data[14]);
	readff  = int(data[15]);
	printf("---------------------------------\n");
	//
	nd=ND;
	ndm=ND-1;
	nd2=ND/2;
	
	//double ch[(nstep+1)*nstep]; // [(nstep+1)*nstep]=[nstep*nstep+nstep], a[i][j]= a[i*n+j] for a[][n]
	//a[z * ySize * xSize + y * xSize + x], a[i][j][k]=a[i*n*m + j*n + k]
	double *ck  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//diffusion potential
	double *ch  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//Concentration data array in tissue
	double *ch2 = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	
	//printf("DELT(0.005)=  "); scanf(" %lf",&delt);	//delt=0.005;
	//printf("ca (0.4)= "); scanf(" %lf",&c0);	//c0=0.4;
	//printf("Mx(1.0) = ");	scanf(" %lf",&Mx);	//Mobility, Mx=0.01;	//Standard value
	//printf("My(1.0) = ");	scanf(" %lf",&My);	//Mobility, My=0.01;	//Standard value
	//printf("Mz(1.0) = ");	scanf(" %lf",&Mz);	//Mobility, Mz=1.0; 	//Standard value

	//al=64.0;		// [nm]
	al=al*1.0e-9;	// [m]

	b1=al/nd;
	//amob_c=1.0;

	time1=0.0; 
	//time1max=30001.0;

	//temp=1000.0;	// [K]

	////L0=2.5e+04/RR/temp;
	//L0=1.0e+03/RR/temp;
	//kapa_c=5.0e-15/b1/b1/RR/temp;
	L0=L00/RR/temp;
	kapa_c=kapa_c0/b1/b1/RR/temp;
	
	//Mx=1.0;  						//mobility
	//My=Mx;  						//mobility
	//Mz=Mx;  						//mobility
	//My=0.0;  						//mobility
	//Mz=0.0;  						//mobility

	//c_flu=c_flu0=0.0;				//Fluctuation coefficient of the field, e.g., 0.0 or s0.1
	c_flu0 = c_flu;

//*************************************************************************
	if(readff==0){
		printf("make initial concentration (random) \n");
		ini_field(ch, ND);
	} else {
		printf("read data.dat file \n");
		datin(ch, ND);
	}
	//gwinsize(INXY,INXY); ginit(2); gsetorg(0,0);//Drawing window display

//**** start **************************************************************
//Nstep = 200;
iout = -1;
start: ;
	printf("time: %f \n", time1);
	if((((int)(time1) % Nstep)==0)) {iout = iout + 1;}
	if((((int)(time1) % Nstep)==0)) {datsave(ch,ND);}
	if((((int)(time1) % Nstep)==0)) {datsave_paraview(ch,ND);}

	//if(time1<1.0e4){c_flu=c_flu0;} else{c_flu=0.0;}

//******[Calculation of diffusion potential********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){

				ip=i+1; im=i-1; jp=j+1; jm=j-1; kp=k+1; km=k-1;
				if(i==ndm) {ip=0;}  if(i==0) {im=ndm;}
				if(j==ndm) {jp=0;}  if(j==0) {jm=ndm;}
				if(k==ndm) {kp=0;}  if(k==0) {km=ndm;}

				//c=ch[i][j][k];
				c=ch[i*ND*ND+j*ND+k];

			 	mu_chem=32.0*L0*c*(1.0-c)*(1.0-2.0*c); 		//chemical potential
			 	//mu_chem=L0*(1.0-2.0*c)+log(c)-log(1.0-c); //chemical potential

				//mu_surf=-2.*kapa_c*(ch[ip][j][k]+ch[im][j][k]
      			//				   +ch[i][jp][k]+ch[i][jm][k]
				//				   +ch[i][j][kp]+ch[i][j][km]-6.0*c);
				mu_surf=-2.*kapa_c*(ch[ip*ND*ND+j*ND+k]+ch[im*ND*ND+j*ND+k]
								   +ch[i*ND*ND+jp*ND+k]+ch[i*ND*ND+jm*ND+k]
      							   +ch[i*ND*ND+j*ND+kp]+ch[i*ND*ND+j*ND+km]
								   -6.0*c);
				
				//ck[i][j][k]=mu_chem+mu_surf;
				ck[i*ND*ND+j*ND+k]= mu_chem+mu_surf;
			}
		}
	}

//******[Time variation of concentration field]**************************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				ip=i+1; im=i-1; jp=j+1; jm=j-1; kp=k+1; km=k-1;
				if(i==ndm) {ip=0;}  if(i==0) {im=ndm;}
				if(j==ndm) {jp=0;}  if(j==0) {jm=ndm;}
				if(k==ndm) {kp=0;}  if(k==0) {km=ndm;}

				//cddtt=Mx*(ck[ip][j][k]+ck[im][j][k]-2.0*ck[i][j][k])
				//     +My*(ck[i][jp][k]+ck[i][jm][k]-2.0*ck[i][j][k])
				//     +Mz*(ck[i][j][kp]+ck[i][j][km]-2.0*ck[i][j][k]);
				cddtt=Mx*(ck[ip*ND*ND+j*ND+k]+ck[im*ND*ND+j*ND+k]-2.0*ck[i*ND*ND+j*ND+k])
					 +My*(ck[i*ND*ND+jp*ND+k]+ck[i*ND*ND+jm*ND+k]-2.0*ck[i*ND*ND+j*ND+k])
					 +Mz*(ck[i*ND*ND+j*ND+kp]+ck[i*ND*ND+j*ND+km]-2.0*ck[i*ND*ND+j*ND+k]);

				//ch2[i][j][k]=ch[i][j][k]+(cddtt+c_flu*(2.0*DRND(1.0)-1.0))*delt;
				ch2[i*ND*ND+j*ND+k] = ch[i*ND*ND+j*ND+k] + (cddtt+c_flu*(2.0*DRND(1.0)-1.0))*delt;
			}
		}
	}

//******[Concentration field balance correction]**********************************
	sumc=0.;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//sumc+=ch2[i][j][k];
				sumc+=ch2[i*ND*ND+j*ND+k];
			}
	   }
	}
    dc0=sumc/nd/nd/nd-c0;

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
	    	for(k=0;k<=ndm;k++){
				//ch[i][j][k]=ch2[i][j][k]-dc0;
				//if(ch[i][j][k]>=1.0){ch[i][j][k]=1.0-1.0e-06;}
				//if(ch[i][j][k]<=0.0){ch[i][j][k]=1.0e-06;}
				ch[i*ND*ND+j*ND+k]=ch2[i*ND*ND+j*ND+k]-dc0;
				if(ch[i*ND*ND+j*ND+k]>=1.0){ch[i*ND*ND+j*ND+k]=1.0-1.0e-06;}
				if(ch[i*ND*ND+j*ND+k]<=0.0){ch[i*ND*ND+j*ND+k]=1.0e-06;}
	    	}
	    }
	}

//******[time increase]*************************************************
	//if(keypress()){return 0;}
	time1=time1+1.;
	if (time1<time1max) {goto start;}

end:;
  return 0;
}

//************[initial concentration wave]*****************************************
void ini_field(double *ch, int ND)
{
	int i, j, k, ii, jj, kk, PN, ipn;
	int x1, y1, z1, r0;
	//int flg;
	int nd=ND, ndm=ND-1, nd2=ND/2; 	//Organization split, organization split-1
	double rnd0, sumc, r; 
 	srand(time(NULL)); // random number initialization
	int nd01, nd02, dd1;

//Initialization
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
				for(k=0;k<=ndm;k++){
					//ch[i][j][k]=0.0;
					ch[i*ND*ND+j*ND+k]=0.0;
			}
		}
	}


// 1: Spinodal decomposition
// 2: Nucleation-growth decomposition
// 3: 1 sphere (volume fraction is calculated here)
// 4: 1 torus (volume fraction is calculated here)

	//flg=1;//Select one of the above.

	switch(flg){

	case 1:	//Spinodal decomposition
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
				for(k=0;k<=ndm;k++){
					//ch[i][j][k]=c0+0.1*(2.*DRND(1)-1.);
					ch[i*ND*ND+j*ND+k]=c0+0.1*(2.*DRND(1)-1.);
			}
		}
	}
	break;


	case 2:	//nucleation
	r0=nd/20;	//put a nucleus of diameter r0
	//PN=50;	//place PN nuclei
	PN=c0/( 4./3.*PI*(1./20.)*(1./20.)*(1./20.) );

	for(ipn=1;ipn<=PN;ipn++){
		x1=nd*DRND(1); y1=nd*DRND(1); z1=nd*DRND(1);
		//while(ch[x1][y1][z1]!=0.0){x1=nd*DRND(1); y1=nd*DRND(1); z1=nd*DRND(1);}
		while(ch[x1*ND*ND+y1*ND+z1]!=0.0){x1=nd*DRND(1); y1=nd*DRND(1); z1=nd*DRND(1);}
		for(i=-r0;i<=(ndm+r0);i++){
			ii=i; if(i<0){ii=nd+i;}  if(i>ndm){ii=i-nd;}
			for(j=-r0;j<=(ndm+r0);j++){
				jj=j; if(j<0){jj=nd+j;}  if(j>ndm){jj=j-nd;}
				for(k=-r0;k<=(ndm+r0);k++){
					kk=k; if(k<0){kk=nd+k;}  if(k>ndm){kk=k-nd;}
					r=sqrt( (double(i-x1))*(double(i-x1))
						   +(double(j-y1))*(double(j-y1))
						   +(double(k-z1))*(double(k-z1)) );
					//if(r<=r0){ ch[ii][jj][kk]=1.0; } 
					if(r<=r0){ ch[ii*ND*ND+jj*ND+kk]=1.0; } 
				}
			}
		}
	}

	sumc=0.;
	//for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ for(k=0;k<=ndm;k++){ sumc+=ch[i][j][k]; } } }
	for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ for(k=0;k<=ndm;k++){ sumc+=ch[i*ND*ND+j*ND+k]; } } }
		//c0=sumc/nd/nd/nd;
	break;


	case 3: //1 ball
	nd01=nd/3;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				dd1=(i-nd2)*(i-nd2)+(j-nd2)*(j-nd2)+(k-nd2)*(k-nd2);
				//if(dd1<(nd01*nd01)){ch[i][j][k]=0.9;} else{ch[i][j][k]=0.2;}
				if(dd1<(nd01*nd01)){ch[i*ND*ND+j*ND+k]=0.9;} else{ch[i*ND*ND+j*ND+k]=0.2;}
			}
		}
	}
	sumc=0.;
	//for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ for(k=0;k<=ndm;k++){ sumc+=ch[i][j][k]; } } }
	for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ for(k=0;k<=ndm;k++){ sumc+=ch[i*ND*ND+j*ND+k]; } } }
		c0=sumc/nd/nd/nd;
	break;


	case 4: //1 torus
	nd01=nd/3; nd02=nd/6;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=nd2-4;k<=nd2+4;k++){
				dd1=(i-nd2)*(i-nd2)+(j-nd2)*(j-nd2);
				//if(dd1<(nd01*nd01)){ch[i][j][k]=0.9;} else{ch[i][j][k]=0.2;}
				//if(dd1<(nd02*nd02)){ch[i][j][k]=0.2;}
				if(dd1<(nd01*nd01)){ch[i*ND*ND+j*ND+k]=0.9;} else{ch[i*ND*ND+j*ND+k]=0.2;}
				if(dd1<(nd02*nd02)){ch[i*ND*ND+j*ND+k]=0.2;}
			}
		}
	}
	sumc=0.;
	//for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ for(k=0;k<=ndm;k++){ sumc+=ch[i][j][k]; } } }
	for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ for(k=0;k<=ndm;k++){ sumc+=ch[i*ND*ND+j*ND+k]; } } }
		c0=sumc/nd/nd/nd;
	break;
	}//switch
}

//************[Saving concentration wave data]************************************
void datsave(double *ch, int ND)
{
	FILE	*stream;
	char	fName[256];
	int 	i, j, k;
	int ndm=ND-1;

	sprintf(fName,"data_%06d.dat",iout);
	//stream = fopen("test3D_lamella.dat", "a");
	stream = fopen(fName, "w");
	fprintf(stream, "%d %d %d \n", ndm, ndm, ndm);
	fprintf(stream, "%e\n", time1);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//fprintf(stream, "%e  ", ch[i][j][k]);
				fprintf(stream, "%e  ", ch[i*ND*ND+j*ND+k]);
			}
		}
	}
	fprintf(stream, "\n");
	fclose(stream);
}

void datsave_paraview(double *ch, int ND)
{
	FILE	*fp;
	char	fName[256];
	int 	i, j, k;
	int ndm=ND-1;
	
	sprintf(fName,"lamella_3D_result%06d.vtk",iout);
	fp = fopen(fName, "w");
	fprintf(fp,"# vtk DataFile Version 3.0 \n");
	fprintf(fp,"output.vtk \n");
	fprintf(fp,"ASCII \n");
	fprintf(fp,"DATASET STRUCTURED_POINTS \n");
	fprintf(fp,"DIMENSIONS %5d %5d %5d \n",(ndm+1),(ndm+1),(ndm+1));
	fprintf(fp,"ORIGIN 0.0 0.0 0.0 \n");
	fprintf(fp,"ASPECT_RATIO 1 1 1 \n");
	fprintf(fp,"POINT_DATA %16d \n",((ndm+1)*(ndm+1)*(ndm+1)));
	fprintf(fp,"SCALARS concentration float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndm;k++){
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				//fprintf(fp,"%10.6f\n", ch[i][j][k]);
				fprintf(fp,"%10.6f\n", ch[i*ND*ND+j*ND+k]);
			}
		}
	}
	fclose(fp);
}
//************ Reading field data *****************************************
void datin(double *ch, int ND)
{
	FILE *datin0;//Stream pointer setting
	int  i, j, k;//integer
	double c00;//Average value of the field
	int nd=ND, ndm=ND-1, nd2=ND/2;
	int rndxm, rndym, rndzm;

	datin0 = fopen("data.dat", "r");//Open the source file

	fscanf(datin0, "%d %d %d", &rndxm, &rndym, &rndzm);
	if (ndm != rndxm){
		printf("data size is mismatch \n");
		printf("Please, change ND= %d in parameters.txt \n", rndxm);
	}
	
	fscanf(datin0, "%lf", &time1);

	c00=0.0;//Initial value of field mean
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				fscanf(datin0, "%lf", &ch[i*ND*ND+j*ND+k]);//Read field data
				c00+=ch[i*ND*ND+j*ND+k];
			}
		}
	}
	c0=c00/nd/nd/nd;//Average value of the field
	printf("c0=  %f  \n", c0);

	fclose(datin0);
}
