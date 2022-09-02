//Read data No.1-5
//Grain number is 1-5
//Consideration of misorientation dependence in interfacial energy density (possible but not done)
//Concentration field is also taken into account (quadratic equation), É¿ phase crystal grain is 1 No.5

#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

#define DRND(x) ((double)(x)/RAND_MAX*rand())//Random number setting

//#define ND 100			//Number of divisions per side of calculation domain in difference calculation
//#define N 6				//Number of crystal orientations to consider + 1

	//int nd=ND, ndm=ND-1;	//Define the number of difference divisions
							// (number of difference blocks) on one side of the computational domain, ND-1
	//int nm=N-1, nmm=N-2;	//Define the number of crystal orientations to consider,
							// N-2 (number of crystal orientations to consider - 1)
	double PI=3.141592;		//pi
	double RR=8.3145;		//gas constant
	double time1;			//calculation count number
	//double ph[N][ND][ND], ph2[N][ND][ND];	//phase field, phase field auxiliary array
	double c0;				//average composition
	//double ch[ND][ND], ch2[ND][ND];	//concentration field, auxiliary array of concentration field
	//double aij[N][N];		//gradient energy factor
	//double wij[N][N];		//the coefficient of the penalty term
	//double tij[N][N];		//Grain boundary mobility
	//double eij[N][N];		//Driving force of grain boundary migration
	//int m00h[N][ND][ND];	//number of directions where p is not 0 at location (i,j) and its surroundings (iÅ}1,jÅ}1)
	//int n00h[ND][ND];		//Number of orientations where p is not 0 at location (i,j) and its surroundings (iÅ}1,jÅ}1)
	int Nstep, iout;

	void ini000(double *ph, double *ch, int ND, int N);	//Initial field setup subroutine
	void datsave(double *ph, double *ch, int ND, int N);//data save subroutine
	void datsave_paraview(double *ph, double *ch, int ND, int N);	//data save subroutine
	void datin(double *ph, int ND, int N);	//data entry subroutine

//******* main program ******************************************
int main(void)
{
	int ND, N;
	int nd, ndm, nm, nmm;
	
	int i, j, k, l, ii, jj, kk, ll, it;	//integer
	int ip, im, jp, jm;	//integer
	int n1, n2, n3;		//integer
	int nalph;
	int n00;			//Number of orientations where p is not 0 at location (i,j) and its surroundings (iÅ}1,jÅ}1)
	//int n000;			//Number of orientations where p is not 0 at position (i,j) (n00>=n000)
	//double cah[ND][ND], cbh[ND][ND];//Concentration field in local equilibrium (parallel tangent rule)
	//double sah[ND][ND], sbh[ND][ND];//Probability of existence of Éø phase and É¿ phase (SÉø and SÉ¿)
	double time1max;	//Maximum calculation count (calculation end count)
	double delt, L, b1;	//Time step, length of one side of calculation area, length of one side of difference block
	double M1, M0;		//Grain boundary mobility
	double W1, W0;		//the coefficient of the penalty term
	double K1, K0;		//gradient energy factor
	double E1, E0;		//Driving force of grain boundary migration
	double temp;		//temperature
	double sum1, sum2, sum3;//Work variables for various sums
	double pddtt;		//Time change rate of phase field

	double gamma, gamma0;//grain boundary energy density
	double delta;		//Grain boundary width (expressed by the number of differential blocks)
	double amobi;		//Grain boundary mobility
	double vm0;			//molar volume

	double A_alph, A_beta, c_alph0, c_beta0;//âªparameters in the theoretical free energy
	double A_alph0, A_beta0;		//parameters in chemical free energy
	double c, sa, sb, Da, Db, D;	//Concentration field, SÉø, SÉ¿, diffusion coefficients (DÉø, DÉ¿), weighted diffusion coefficient D
	double gii, gjj, cii, cjj;		//chemical free energy, local equilibrium concentration
	double dgdc;					//Concentration derivative of chemical free energy
	double cddtt;					//concentration change rate
	double dc0;						//Working variables for density field correction
	double dev1_a, dev1_b, dev2_a, dev2_b;//working variables in the differential calculation in the diffusion equation


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
	nalph   = data[15];
	A_alph0 = data[16];
	c_alph0 = data[17];
	A_beta0 = data[18];
	c_beta0 = data[19];
	Da      = data[20];
	Db      = data[21];
	D       = data[22];
	printf("---------------------------------\n");
	//
	nd=ND;
	ndm=ND-1;
	nm=N-1;
	nmm=N-2;
	//
	double *ph   = (double *)malloc(sizeof(double)*( N*ND*ND + ND*ND + ND ));	//phase field
	double *ph2  = (double *)malloc(sizeof(double)*( N*ND*ND + ND*ND + ND ));	//Phase field auxiliary array
	double *ch   = (double *)malloc(sizeof(double)*( ND*ND + ND ));
	double *ch2  = (double *)malloc(sizeof(double)*( ND*ND + ND ));
	double *aij  = (double *)malloc(sizeof(double)*( ND*ND + ND ));	//gradient energy factor
	double *wij  = (double *)malloc(sizeof(double)*( ND*ND + ND ));	//the coefficient of the penalty term
	double *tij  = (double *)malloc(sizeof(double)*( ND*ND + ND ));	//Grain boundary mobility
	double *eij  = (double *)malloc(sizeof(double)*( ND*ND + ND ));	//Driving force of grain boundary migration
	double *m00h = (double *)malloc(sizeof(double)*( N*ND*ND + ND*ND + ND ));	//number of orientations where p is not 0 at location (i,j) and its surroundings (iÅ}1,jÅ}1)
	double *n00h = (double *)malloc(sizeof(double)*( ND*ND + ND ));	//Number of orientations where p is not 0 at location (i,j) and its surroundings (iÅ}1,jÅ}1)
	//
	double *cah  = (double *)malloc(sizeof(double)*( ND*ND + ND ));	//Concentration field in local equilibrium (parallel tangent rule)
	double *cbh  = (double *)malloc(sizeof(double)*( ND*ND + ND ));	//Concentration field in local equilibrium (parallel tangent rule)
	double *sah  = (double *)malloc(sizeof(double)*( ND*ND + ND ));	//Probability of existence of Éø phase and É¿ phase (SÉø and SÉ¿)
	double *sbh  = (double *)malloc(sizeof(double)*( ND*ND + ND ));	//Probability of existence of Éø phase and É¿ phase (SÉø and SÉ¿)
	//
	//printf("delt(0.2)=  "); scanf(" %lf",&delt);//Enter time step	//	delt=0.2;
	//printf("c0(0.4)=  "); scanf(" %lf",&c0);//Enter Average Composition	//c0=0.2;

	//temp=1000.0; 				//temperature (K)
	//L=2000.0;					//Side length of calculation area (nm)
	b1=L/(double)ND*1.0e-9; 	//Length of one side of difference block (m)

	//vm0=7.0e-6;				//molar volume
	//gamma=0.5*vm0/RR/temp/b1;	//Dimensionless grain boundary energy density (0.5J/m^2)
	gamma=gamma0*vm0/RR/temp/b1;//Dimensionless grain boundary energy density (0.5J/m^2)
	//delta=5.0;				//Grain boundary width (expressed by the number of differential blocks)

	//K1=8.0*delta*gamma/PI/PI;	//gradient energy coefficient [equation (3.23)]
	K1=K0*delta*gamma/PI/PI;	//gradient energy coefficient [equation (3.23)]
	//W1=4.0*gamma/delta;		//the coefficient of the penalty term [equation (3.23)]
	W1=W0*gamma/delta;			//the coefficient of the penalty term [equation (3.23)]
	//amobi=1.0;
	//M1=amobi*PI*PI/(8.0*delta);	//Grain boundary mobility [equation (3.23)]
	M1=amobi*PI*PI/(M0*delta);	//Grain boundary mobility [equation (3.23)]
	//E1=50.0/RR/temp;			//Driving force of grain boundary migration
	E1=E0/RR/temp;				//Driving force of grain boundary migration

	time1=0.0;
	//time1max=10000001.;		//Initial and maximum calculation counts

	//nalph=4;					//Number of Éø-phase crystal orientations (number of Éø-phase grains considered)
	//A_alph=1.0e+02/RR/temp;	//parameters in chemical free energy
	A_alph=A_alph0/RR/temp;		//parameters in chemical free energy
	//c_alph0=0.1;
	//A_beta=1.0e+02/RR/temp;
	A_beta=A_beta0/RR/temp;
	//c_beta0=0.9;
	//Da=Db=D=1.0;				//Assuming equal diffusion coefficients for each phase

//*** Setting the array (K,W,M,E) of Eq.(4.36)-Eq.(4.39) **************************
	for(i=1;i<=nm;i++){
		for(j=i+1;j<=nm;j++){
			//wij[i][j]=wij[j][i]=W1;
			//aij[i][j]=aij[j][i]=K1;
			//tij[i][j]=tij[j][i]=M1;
			//eij[i][j]=eij[j][i]=0.0;
			wij[i*ND+j]=wij[j*ND+i]=W1;
			aij[i*ND+j]=aij[j*ND+i]=K1;
			tij[i*ND+j]=tij[j*ND+i]=M1;
			eij[i*ND+j]=eij[j*ND+i]=0.0;
		}
		//wij[i][i]=0.0;  aij[i][i]=0.0;  tij[i][i]=0.0;  eij[i][i]=0.0;
		wij[i*ND+i]=0.0; aij[i*ND+i]=0.0; tij[i*ND+i]=0.0; eij[i*ND+i]=0.0;
	}

//*** Initial Field Settings and Drawing Window Display *****************************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=1;k<=nm;k++){
				//ph[k][i][j]=0.0;  ph2[k][i][j]=0.0;//Phase field and auxiliary array initialization
				 ph[k*ND*ND+i*ND+j]=0.0;	//Phase field initialization
				ph2[k*ND*ND+i*ND+j]=0.0;	//Auxiliary array initialization
			}
			//ch[i][j]=0.0;  ch2[i][j]=0.0;//Initialization of the concentration field and auxiliary arrays
			 ch[i*ND+j]=0.0;	//Initialization of the concentration field
			ch2[i*ND+j]=0.0;	//Auxiliary array initialization
		}
	}

	ini000(ph, ch, ND, N);//Initial field setting
	//datin(ph, ND, N);//Initial field input

//**** Simulation start ******************************
//Nstep = 200;
iout = -1;
start: ;

	if((((int)(time1) % Nstep)==0)) {datsave(ph, ch, ND, N);}//Save the field every fixed repetition count
	if((((int)(time1) % Nstep)==0)) {datsave_paraview(ph, ch, ND, N);}//Save the field every fixed repetition count
	//if(time1==200.) {datsave();}//save the field for a specific time

//****** Calculation of local equilibrium composition (cÉø and cÉ¿) [equation (4.47)] ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//c=ch[i][j];	//concentration field
			c=ch[i*ND+j];	//concentration field
			//sum1=0.0; for(ii=1;ii<=nalph;ii++){ sum1+=ph[ii][i][j]; }
			sum1=0.0; for(ii=1;ii<=nalph;ii++){ sum1+=ph[ii*ND*ND+i*ND+j]; }
			//sa=sah[i][j]=sum1;	//Calculation of SÉø [equation (4.46)]
			sa=sah[i*ND+j]=sum1;	//Calculation of SÉø [equation (4.46)]
			//sum1=0.0; for(ii=nalph+1;ii<=nm;ii++){ sum1+=ph[ii][i][j]; }
			sum1=0.0; for(ii=nalph+1;ii<=nm;ii++){ sum1+=ph[ii*ND*ND+i*ND+j]; }
			//sb=sbh[i][j]=sum1;	//Calculation of SÉ¿ [equation (4.46)]
			sb=sbh[i*ND+j]=sum1;	//Calculation of SÉ¿ [equation (4.46)]

			//Calculation of local equilibrium composition [equation (4.47)]
			//cah[i][j]=(A_beta*c+(A_alph*c_alph0-A_beta*c_beta0)*sb)/(A_beta*sa+A_alph*sb);
			cah[i*ND+j]=(A_beta*c+(A_alph*c_alph0-A_beta*c_beta0)*sb)/(A_beta*sa+A_alph*sb);
			//cbh[i][j]=(A_alph*c+(A_beta*c_beta0-A_alph*c_alph0)*sa)/(A_beta*sa+A_alph*sb);
			cbh[i*ND+j]=(A_alph*c+(A_beta*c_beta0-A_alph*c_alph0)*sa)/(A_beta*sa+A_alph*sb);
			//if(cah[i][j]>=1.0){cah[i][j]=1.0;}  if(cah[i][j]<=0.0){cah[i][j]=0.0;}//Domain correction of the concentration field
			if(cah[i*ND+j]>=1.0){cah[i*ND+j]=1.0;}  if(cah[i*ND+j]<=0.0){cah[i*ND+j]=0.0;}//Domain correction of the concentration field
			//if(cbh[i][j]>=1.0){cbh[i][j]=1.0;}  if(cbh[i][j]<=0.0){cbh[i][j]=0.0;}
			if(cbh[i*ND+j]>=1.0){cbh[i*ND+j]=1.0;}  if(cbh[i*ND+j]<=0.0){cbh[i*ND+j]=0.0;}
		}
	}

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
				//  ( (ph[ii][i][j]==0.0)&&(ph[ii][ip][j]>0.0)||
				//  	(ph[ii][im][j]>0.0)||
				//  	(ph[ii][i][jp]>0.0)||
				//  	(ph[ii][i][jm]>0.0)   ) ){
				//  		n00++; m00h[n00][i][j]=ii;
				//		//printf("%d  ", n00);
				//}
				if( (ph[ii*ND*ND+i*ND+j]>0.0)||
				  ( (ph[ii*ND*ND+i*ND+j]==0.0)&&(ph[ii*ND*ND+ip*ND+j]>0.0)||
				  	(ph[ii*ND*ND+im*ND+j]>0.0)||
				  	(ph[ii*ND*ND+i*ND+jp]>0.0)||
				  	(ph[ii*ND*ND+i*ND+jm]>0.0)   ) ){
				  		n00++; m00h[n00*ND*ND+i*ND+j]=ii;
						//printf("%d  ", n00);
				}
			}
			//n00h[i][j]=n00;
			n00h[i*ND+j]=n00;
//--------------------------------------------------------------------------
		}
	}


//***** Calculation of Phase Field Evolution Equation **********************************
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
						//		 +(wij[ii][kk]-wij[jj][kk])*ph[kk][i][j];//[part of formula (4.41)]
						sum1+=0.5*(aij[ii*ND+kk]-aij[jj*ND+kk])*(ph[kk*ND*ND+ip*ND+j]+ph[kk*ND*ND+im*ND+j]
																+ph[kk*ND*ND+i*ND+jp]+ph[kk*ND*ND+i*ND+jm]-4.0*ph[kk*ND*ND+i*ND+j])
								 +(wij[ii*ND+kk]-wij[jj*ND+kk])*(ph[kk*ND*ND+i*ND+j]);//[part of formula (4.41)]
					}

					//if(ii<=nalph){gii=A_alph*(cah[i][j]-c_alph0)*(cah[i][j]-c_alph0); cii=cah[i][j]; }
					if(ii<=nalph){gii=A_alph*(cah[i*ND+j]-c_alph0)*(cah[i*ND+j]-c_alph0); cii=cah[i*ND+j]; }
					//else{gii=A_beta*(cbh[i][j]-c_beta0)*(cbh[i][j]-c_beta0); cii=cbh[i][j]; }//formula (4.43)
					else{gii=A_beta*(cbh[i*ND+j]-c_beta0)*(cbh[i*ND+j]-c_beta0); cii=cbh[i*ND+j]; }//formula (4.43)
					//if(jj<=nalph){gjj=A_alph*(cah[i][j]-c_alph0)*(cah[i][j]-c_alph0); cjj=cah[i][j]; }
					if(jj<=nalph){gjj=A_alph*(cah[i*ND+j]-c_alph0)*(cah[i*ND+j]-c_alph0); cjj=cah[i*ND+j]; }
					//else{gjj=A_beta*(cbh[i][j]-c_beta0)*(cbh[i][j]-c_beta0); cjj=cbh[i][j]; }//formula (4.43)
					else{gjj=A_beta*(cbh[i*ND+j]-c_beta0)*(cbh[i*ND+j]-c_beta0); cjj=cbh[i*ND+j]; }//formula (4.43)
						//dgdc=2.0*A_alph*(cah[i][j]-c_alph0);//Calculation of concentration derivative of chemical free energy [equation (4.48)]
						dgdc=2.0*A_alph*(cah[i*ND+j]-c_alph0);//Calculation of concentration derivative of chemical free energy [equation (4.48)]
						//dgdc=2.0*A_beta*(cbh[i][j]-c_beta0);
						
						//pddtt+=-2.0*tij[ii][jj]/double(n00h[i][j])*(sum1+gii-gjj-(cii-cjj)*dgdc );
						pddtt+=-2.0*tij[ii*ND+jj]/double(n00h[i*ND+j])*(sum1+gii-gjj-(cii-cjj)*dgdc );
						//hase-field evolution equation [equation (4.41)]
				}
				//ph2[ii][i][j]=ph[ii][i][j]+pddtt*delt;//Phase field time evolution (explicit method)
				//if(ph2[ii][i][j]>=1.0){ph2[ii][i][j]=1.0;}//Phase field domain correction
				//if(ph2[ii][i][j]<=0.0){ph2[ii][i][j]=0.0;}
				ph2[ii*ND*ND+i*ND+j]=ph[ii*ND*ND+i*ND+j]+pddtt*delt;//Phase field time evolution (explicit method)
				if(ph2[ii*ND*ND+i*ND+j]>=1.0){ph2[ii*ND*ND+i*ND+j]=1.0;}//Phase field domain correction
				if(ph2[ii*ND*ND+i*ND+j]<=0.0){ph2[ii*ND*ND+i*ND+j]=0.0;}
			}
		}//j
	}//i

//***** Calculation of time evolution of concentration field (diffusion equation) **********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ip=i+1;  im=i-1;  jp=j+1;  jm=j-1;
			if(i==ndm){ip=0;}  if(i==0){im=ndm;}//periodic boundary conditions
			if(j==ndm){jp=0;}  if(j==0){jm=ndm;}

			//Differential calculation in the diffusion equation
			//dev1_a=0.25*( (sah[ip][j]-sah[im][j])*(cah[ip][j]-cah[im][j])
			//			 +(sah[i][jp]-sah[i][jm])*(cah[i][jp]-cah[i][jm]) );
			dev1_a=0.25*( (sah[ip*ND+j]-sah[im*ND+j])*(cah[ip*ND+j]-cah[im*ND+j])
						 +(sah[i*ND+jp]-sah[i*ND+jm])*(cah[i*ND+jp]-cah[i*ND+jm]) );
			//dev1_b=0.25*( (sbh[ip][j]-sbh[im][j])*(cbh[ip][j]-cbh[im][j])
			//			 +(sbh[i][jp]-sbh[i][jm])*(cbh[i][jp]-cbh[i][jm]) );
			dev1_b=0.25*( (sbh[ip*ND+j]-sbh[im*ND+j])*(cbh[ip*ND+j]-cbh[im*ND+j])
						 +(sbh[i*ND+jp]-sbh[i*ND+jm])*(cbh[i*ND+jp]-cbh[i*ND+jm]) );
			//dev2_a=sah[i][j]*(cah[ip][j]+cah[im][j]+cah[i][jp]+cah[i][jm]-4.0*cah[i][j]);
			dev2_a=sah[i*ND+j]*(cah[ip*ND+j]+cah[im*ND+j]+cah[i*ND+jp]+cah[i*ND+jm]-4.0*cah[i*ND+j]);
			//dev2_b=sbh[i][j]*(cbh[ip][j]+cbh[im][j]+cbh[i][jp]+cbh[i][jm]-4.0*cbh[i][j]);
			dev2_b=sbh[i*ND+j]*(cbh[ip*ND+j]+cbh[im*ND+j]+cbh[i*ND+jp]+cbh[i*ND+jm]-4.0*cbh[i*ND+j]);

			cddtt=Da*(dev1_a+dev2_a)+Db*(dev1_b+dev2_b);	//Diffusion equation [equation (4.42)]
			//ch2[i][j]=ch[i][j]+cddtt*delt;	//Time evolution of concentration field (explicit method)
			//ch2[i][j]=ch[i][j]+cddtt*delt+(2.*DRND(1.)-1.)*0.001;	//Time evolution of concentration field (explicit method)
			ch2[i*ND+j]=ch[i*ND+j]+cddtt*delt+(2.*DRND(1.)-1.)*0.001;	//Time evolution of concentration field (explicit method)
		}
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=1;k<=nm;k++){
				//ph[k][i][j]=ph2[k][i][j];//Move auxiliary array to main array (phase field)
				ph[k*ND*ND+i*ND+j]=ph2[k*ND*ND+i*ND+j];//Move auxiliary array to main array (phase field)
			}
			//ch[i][j]=ch2[i][j];//Move auxiliary array to main array (concentration field)
			ch[i*ND+j]=ch2[i*ND+j];//Move auxiliary array to main array (concentration field)
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

//*** Concentration field balance correction *************************************************************
	//sum1=0.0; for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ sum1+=ch[i][j]; } }
	sum1=0.0; for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ sum1+=ch[i*ND+j]; } }
	dc0=sum1/nd/nd-c0;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//ch[i][j]=ch[i][j]-dc0;
			//if(ch[i][j]>1.0){ch[i][j]=1.0;}   if(ch[i][j]<0.0){ch[i][j]=0.0;}
			ch[i*ND+j]=ch[i*ND+j]-dc0;
			if(ch[i*ND+j]>1.0){ch[i*ND+j]=1.0;}
			if(ch[i*ND+j]<0.0){ch[i*ND+j]=0.0;}
	  }
	}

//*********************************************************************
	//if(keypress()){return 0;}//Waiting for key
	time1=time1+1.;
	
	if(time1<time1max){goto start;}//Determining if the maximum count has been reached
	printf("Finished \n");

end:;
	std::exit(0);
	//return 0;
}

//************ Initial field (phase field and concentration field) setting subroutine *************
void ini000(double *ph, double *ch, int ND, int N)
{
	int i, j, k, l, it;		//integer
	int ii, jj, kk;			//integer
	int ip, im, jp, jm;		//integer
	int igx, igy, ixmin, ixmax, iymin, iymax;	//Setting the Difference Block Coordinate System
	double x, y, xmin, xmax, ymin, ymax;		//Setting the normalized coordinate system
	double sum1, sum2, t, r0, phi, r;			//work variable
	//srand(time(NULL)); //random number initialization
	int nd=ND, ndm=ND-1, nm=N-1, nmm=N-2;

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//ch[i][j]=0.5;//Setting the Éø-phase concentration field (supersaturated solid solution)
			ch[i*ND+j]=0.5;//Setting the Éø-phase concentration field (supersaturated solid solution)
		}
	}

	datin(ph, ND, N);//Read crystal structure data of Éø phase

	xmin=-1.0; xmax=1.0; ymin=-1.0;  ymax=1.0;	//normalized coordinate system
	ixmin=0; ixmax=ND; iymin=0;  iymax=ND;		//difference block coordinate system

	r0=2.0;//É¿ phase size

//Below, the nuclei of the É¿ phase are placed at the triple points (6 points) of the grain boundaries of the Éø phase.
	x=0.5;           igx=(ixmax-ixmin)/(xmax-xmin)*(x-xmin)+ixmin;
	y=0.5/sqrt(3.);  igy=(iymax-iymin)/(ymax-ymin)*(y-ymin)+iymin;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			r=sqrt( (double(i-igx))*(double(i-igx))+(double(j-igy))*(double(j-igy)) );
			//if(r<=r0){ ph[nm][i][j]=1.0; 	for(k=1;j<=nm-1;k++){ph[k][i][j]=0.0;} ch[i][j]=0.9; }
			if(r<=r0){ ph[nm*ND*ND+i*ND+j]=1.0;
				for(k=1;j<=nm-1;k++){ph[k*ND*ND+i*ND+j]=0.0;} 
				ch[i*ND+j]=0.9;
			}
		}
	}

	x=-0.5;          igx=(ixmax-ixmin)/(xmax-xmin)*(x-xmin)+ixmin;
	y=0.5/sqrt(3.);  igy=(iymax-iymin)/(ymax-ymin)*(y-ymin)+iymin;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			r=sqrt( (double(i-igx))*(double(i-igx))+(double(j-igy))*(double(j-igy)) );
			//if(r<=r0){ ph[nm][i][j]=1.0; 	for(k=1;j<=nm-1;k++){ph[k][i][j]=0.0;} ch[i][j]=0.9;  }
			if(r<=r0){
				ph[nm*ND*ND+i*ND+j]=1.0;
				for(k=1;j<=nm-1;k++){ph[k*ND*ND+i*ND+j]=0.0;}
				ch[i*ND+j]=0.9;
			}
		}
	}

	x=0.5;           igx=(ixmax-ixmin)/(xmax-xmin)*(x-xmin)+ixmin;
	y=-0.5/sqrt(3.); igy=(iymax-iymin)/(ymax-ymin)*(y-ymin)+iymin;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			r=sqrt( (double(i-igx))*(double(i-igx))+(double(j-igy))*(double(j-igy)) );
			//if(r<=r0){ ph[nm][i][j]=1.0; 	for(k=1;j<=nm-1;k++){ph[k][i][j]=0.0;} ch[i][j]=0.9;  }
			if(r<=r0){
				ph[nm*ND*ND+i*ND+j]=1.0;
				for(k=1;j<=nm-1;k++){ph[k*ND*ND+i*ND+j]=0.0;}
				ch[i*ND+j]=0.9;
			}
		}
	}

	x=-0.5;           igx=(ixmax-ixmin)/(xmax-xmin)*(x-xmin)+ixmin;
	y=-0.5/sqrt(3.);  igy=(iymax-iymin)/(ymax-ymin)*(y-ymin)+iymin;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			r=sqrt( (double(i-igx))*(double(i-igx))+(double(j-igy))*(double(j-igy)) );
			//if(r<=r0){ ph[nm][i][j]=1.0; 	for(k=1;j<=nm-1;k++){ph[k][i][j]=0.0;} ch[i][j]=0.9;  }
			if(r<=r0){
				ph[nm*ND*ND+i*ND+j]=1.0;
				for(k=1;j<=nm-1;k++){ph[k*ND*ND+i*ND+j]=0.0;}
				ch[i*ND+j]=0.9;
			}
		}
	}

	x=0.0;          igx=(ixmax-ixmin)/(xmax-xmin)*(x-xmin)+ixmin;
	y=1./sqrt(3.);  igy=(iymax-iymin)/(ymax-ymin)*(y-ymin)+iymin;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			r=sqrt( (double(i-igx))*(double(i-igx))+(double(j-igy))*(double(j-igy)) );
			//if(r<=r0){ ph[nm][i][j]=1.0; 	for(k=1;j<=nm-1;k++){ph[k][i][j]=0.0;} ch[i][j]=0.9;  }
			if(r<=r0){
				ph[nm*ND*ND+i*ND+j]=1.0;
				for(k=1;j<=nm-1;k++){ph[k*ND*ND+i*ND+j]=0.0;}
				ch[i*ND+j]=0.9;
			}
		}
	}

	x=0.0;           igx=(ixmax-ixmin)/(xmax-xmin)*(x-xmin)+ixmin;
	y=-1./sqrt(3.);  igy=(iymax-iymin)/(ymax-ymin)*(y-ymin)+iymin;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			r=sqrt( (double(i-igx))*(double(i-igx))+(double(j-igy))*(double(j-igy)) );
			//if(r<=r0){ ph[nm][i][j]=1.0; 	for(k=1;j<=nm-1;k++){ph[k][i][j]=0.0;} ch[i][j]=0.9;  }
			if(r<=r0){
				ph[nm*ND*ND+i*ND+j]=1.0;
				for(k=1;j<=nm-1;k++){ph[k*ND*ND+i*ND+j]=0.0;}
				ch[i*ND+j]=0.9;
			}
		}
	}

//--- Phase field normalization correction and average composition calculation ------------------------
	sum2=0.0;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//sum1=0.0; for(k=1;k<=nm;k++){ sum1+=ph[k][i][j]; }
			//for(k=1;k<=nm;k++){ ph[k][i][j]=ph[k][i][j]/sum1; }
			//sum2+=ch[i][j];
			sum1=0.0; for(k=1;k<=nm;k++){ sum1+=ph[k*ND*ND+i*ND+j]; }
			for(k=1;k<=nm;k++){ ph[k*ND*ND+i*ND+j]=ph[k*ND*ND+i*ND+j]/sum1; }
			sum2+=ch[i*ND+j];
		}
	}
	c0=sum2/nd/nd;

}

//************ data save subroutine *******************************
void datsave(double *ph, double *ch, int ND, int N)
{
	FILE		*stream;//Stream pointer setting
	int 		i, j, k;//integer
	double 	col;
	int nd=ND, ndm=ND-1, nm=N-1, nmm=N-2;

	stream = fopen("test.dat", "a");//Open the file to write to in append mode
	fprintf(stream, "%e  \n", time1);//Saving calculation counts
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=nm;k++){
				//fprintf(stream, "%e   ", ph[k][i][j]);//Save phase field
				fprintf(stream, "%e   ", ph[k*ND*ND+i*ND+j]);//Save phase field
			}
			//fprintf(stream, "%e   ", ch[i][j]);//Concentration field storage
			fprintf(stream, "%e   ", ch[i*ND+j]);//Concentration field storage
		}
	}
	fprintf(stream, "\n");//writing a newline
	fclose(stream);//close file
}

void datsave_paraview(double *ph, double *ch, int ND, int N)
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
		sprintf(fName,"mpf_N%03d_result%06d.vtk",k,iout);
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
		sprintf(fName,"mpf_N%03d_result%06d.vtk",k,iout);
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
		fprintf(fp,"SCALARS concentration float \n");
		fprintf(fp,"LOOKUP_TABLE default \n");
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				//fprintf(stream, "%e   ", ch[i][j]);//Concentration field storage
				//fprintf(fp,"%10.6f\n", ch[i][j]);
				fprintf(fp,"%10.6f\n", ch[i*ND+j]);
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
	double dami1;//work variable
	int nd=ND, ndm=ND-1, nm=N-1, nmm=N-2;

	datin0 = fopen("test.dat", "r");//open source file
	fscanf(datin0, "%lf  ", &dami1);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=1;k<=nm;k++){
				//fscanf(datin0, "%lf  ", &ph[k][i][j]);//Load phase field
				fscanf(datin0, "%lf  ", &ph[k*ND*ND+i*ND+j]);//Load phase field
			}
		}
	}
	fclose(datin0);//close file

}

