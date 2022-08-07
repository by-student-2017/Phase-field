//*** [program (FePt2D_001_v2.java)] ************************************************
//*** [import] ****************************
import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.util.*; // for parameters.txt

public class FePt2D_001_v2 extends Frame{
//*** [Global variables] *****************************************************************************************
	static int ND;					//Number of divisions on one side of the organization
	static int IG;					//2^IG=ND
	static int nd;					//Number of divisions on one side of the calculation area
	static int ndm;					//Number of divisions on one side of the calculation area-1
	static int nd2;					//Number of divisions on one side of the calculation area/2
	static int ig;					//2^ig=ND
	static int width;				//Overall width of the window
	static int height;				//Overall height of the window
	static int xwidth;				//Drawing area width
	static int yheight;				//Drawing area height
	static int insetx;				//Window frame width (left and right and bottom)
	static int insety;				//Window frame width (top)
	static double PI=3.141592;		//pi
	static double RR=8.3145;		//Gas constant
	static double [][] s1h;			//Order parameter data s1 array in the organization
	static double [][] s2h;			//Order parameter data s2 array in the organization
	static Graphics g;				//Free energy curve screen graphics object
	static double time1;			//Calculation time (count number)
	static double temp; 			//Temperature [K]

	//For the program of Fast Fourier Transform, refer to Reference (9).
	static double qs;				//Distinguishing between Fourier transform (qs: -1) and inverse Fourier transform (qs: 1)
	static double [][] xr;			//Array used for the real part of the Fourier transform
	static double [][] xi;			//Array used for the imaginary part of the Fourier transform
	static double [] xrf;			//Array used for the real part of the Fourier transform
	static double [] xif;			//Array used for the imaginary part of the Fourier transform
	static double [] s;				//sin table
	static double [] c;				//cos table
	static int [] ik;				//Array of bit inversion operations

//*** [constructor] ****************************
	public FePt2D_001_v2(){
		xwidth=400; yheight=400; 			//Horizontal and vertical length of the drawing screen (in pixels)
		insetx=4; insety=30;				//Border length of drawing screen
		width=xwidth+insetx*2;  			//Horizontal length of the entire drawing window
		height=yheight+insetx+insety; 		//Vertical length of the entire drawing window
		setSize(width, height);				//Set of drawing windows
		setBackground(Color.white); 		//Set the color of the drawing part of the drawing window to white
		setVisible(true); 					//Make the drawing window visible
		addWindowListener(new WindowAdapter(){ 
			public void windowClosing(WindowEvent e){ System.exit(0); }
											//Operation when closing Window (Setting of x in the upper right of Window)
		});
	}

//*** [main program] *******************************************************************************************
	public static void main(String[] args) throws Exception{//Exception handling is not performed

		FePt2D_001_v2 prog=new FePt2D_001_v2();		//Generate an instance prog of FePt2D_001_v2

		double [][] ec11;							//Constrained strain fluctuation amount
		double [][] ec22;							//Constrained strain fluctuation amount
		double [][] ep11h0;							//Metamorphosis distortion
		double [][] ep22h0;							//Metamorphosis distortion
		double [][] ep11qrh0;						//Fourier transform of constrained strain fluctuation (real part)
		double [][] ep11qih0;						//Fourier transform of constrained strain fluctuation amount (imaginary part)
		double [][] ep22qrh0;						//Fourier transform of constrained strain fluctuation (real part)
		double [][] ep22qih0;						//Fourier transform of constrained strain fluctuation amount (imaginary part)
		double [][] eta_s1 = new double[4][4];		//Variant 1 transformation strain
		double [][] eta_s2 = new double[4][4];		//Variant 2 transformation strain
		double [][] s1k_su;							//Gradient potential
		double [][] s2k_su;							//Gradient potential

		double s1, s2;								//Martensite phase field
		double s1k_chem, s1k_str;					//Chemical potential, elastic potential
		double s2k_chem, s2k_str;					//Chemical potential, elastic potential
		double c11, c12, c44, lam0, mu0, nu0;		//Elastic constant
		double c11c, c12c, c44c;
		double eta_s1_11, eta_s1_22, eta_s1_33;		//Variant 1 transformation strain
		double eta_s1_12, eta_s1_13, eta_s1_23;		//Variant 1 transformation strain
		double el_fac, el_fac_c;					//Elastic constant normalization variable
		double ep11T, ep22T;						//Sum of individual strain components
		double ep11_0, ep22_0;						//Mean value of transformation strain in tissue
		double ep11_a, ep22_a, ep12_a, ep21_a;		//Distortion due to external force
		double sig11_a, sig22_a;					//External force
		double Z11ep, Z12ep, Z21ep, Z22ep;			//Coefficients used during Fourier transform
		double sum11, sum22;
		double s1ddtt, s2ddtt;						//Amount of time change in phase field
		double delt;								//Time increments

		int i, j, k, l, ii=0, jj=0;					//integer
		int p, q, m, n;								//integer
		int ip, im, jp, jm;							//integer
		int Nstep;									//integer
		double al, temp;							//Calculation area, temperature
		double time1max;							//Maximum time (used to stop calculation)
		double b1, vm0, atom_n;						//Normalized length, molar volume, number of atoms in a unit cell
		double smob;								//Relaxation coefficient of martensitic transformation dynamics
		double nx, ny, nxx, nyy, alnn;				//Basic vector of reciprocal space, its square, norm

		double AA0, AA1, AA2, AA3, AA4, AA5, AA6;	//Chemical driving force constant
		double a1_c, b1_c, c1_c;					//Lattice constant
		double kappa_s1, kappa_s2;					//Gradient energy constant
		double ds_fac;								//The magnitude of the fluctuation of the phase field

		double [] data = new double[40];
		
//---- Various parameter settings ----------------------------------------------------
		System.out.printf("------------------------------------- \n");
		String fileName = "parameters.txt";
		Scanner scan = new Scanner(new File(fileName));
		i = 0;
		while(scan.hasNextLine()){
			String line = scan.nextLine();
			String[] str = line.split(" ");
			//System.out.println(str[1]);
			data[i] = Double.parseDouble(str[1]);
			System.out.printf("%d,%s %f \n",i, str[0],data[i]);
			i = i + 1;
		}
		System.out.printf("------------------------------------- \n");
		ND       = (int)data[0];
		delt     = data[1];
		temp     = data[2];
		al       = data[3];
		time1max = data[4];
		smob     = data[5];
		ds_fac   = data[6];
		AA0      = data[7];
		AA1      = data[8];
		AA2      = data[9];
		AA3      = data[10];
		AA4      = data[11];
		AA5      = data[12];
		AA6      = data[13];
		kappa_s1 = data[14];
		kappa_s2 = data[15];
		vm0      = data[16];
		eta_s1_11= data[17];
		eta_s1_22= data[18];
		eta_s1_33= data[19];
		eta_s1_12= data[20];
		eta_s1_13= data[21];
		eta_s1_23= data[22];
		el_fac_c = data[23];
		c11c     = data[24];
		c12c     = data[25];
		c44c     = data[26];
		sig22_a  = data[27];
		ep11_a   = data[28];
		ep22_a   = data[29];
		ep12_a   = data[30];
		ep21_a   = data[31];
		//
		IG=(int) (Math.log((double) ND) / Math.log(2.0));
		nd=ND;
		ndm=ND-1;	//Number of divisions on one side of the calculation area-1
		nd2=ND/2;	//Number of divisions on one side of the calculation area/2
		ig=IG;		//2^ig=ND
		//
		s1h = new double[ND][ND];	//Order parameter data s1 array in the organization
		s2h = new double[ND][ND];	//Order parameter data s2 array in the organization
		//
		xr  = new double[ND][ND];	//Array used for the real part of the Fourier transform
		xi  = new double[ND][ND];	//Array used for the imaginary part of the Fourier transform
		xrf = new double[ND];		//Array used for the real part of the Fourier transform
		xif = new double[ND];		//Array used for the imaginary part of the Fourier transform
		s   = new double[ND];		//sin table
		c   = new double[ND];		//cos table
		ik  = new int[ND];			//Array of bit inversion operations
		//
		ec11     = new double[ND][ND];	//Constrained strain fluctuation amount
		ec22     = new double[ND][ND];	//Constrained strain fluctuation amount
		ep11h0   = new double[ND][ND];	//Metamorphosis distortion
		ep22h0   = new double[ND][ND];	//Metamorphosis distortion
		ep11qrh0 = new double[ND][ND];	//Fourier transform of constrained strain fluctuation (real part)
		ep11qih0 = new double[ND][ND];	//Fourier transform of constrained strain fluctuation amount (imaginary part)
		ep22qrh0 = new double[ND][ND];	//Fourier transform of constrained strain fluctuation (real part)
		ep22qih0 = new double[ND][ND];	//Fourier transform of constrained strain fluctuation amount (imaginary part)
		//
		s1k_su   = new double[ND][ND];	//Gradient potential
		s2k_su   = new double[ND][ND];	//Gradient potential
		//
		//delt=0.05;					//Enter time

		//temp=873.0;					//Temperature [K]
		//al=200.0*1.0E-09;			//Computation area [m]
		al=al*1.0E-09;			//Computation area [m]
		b1=al/nd;					//Difference block length [m]

		time1=-10.0;				//Initial setting time
		//time1max=1.0+1.0e+07;		//Maximum calculation time

		//smob=1.0;					//Relaxation coefficient of martensitic transformation dynamics
		//ds_fac=0.1;					//Fluctuation coefficient of phase field

		//AA0=3.82027e+03;			//Chemical driving force of martensitic transformation [J/mol]
		AA0=AA0/RR/temp;			//Dimensionless
		//AA1=0.1; 					//Chemical driving force constant
		//AA2=-4.0*AA1-12.0; 
		//AA3=3.0*AA1+12.0;	
		//AA4=AA5=AA6=4.0;

		//kappa_s1=1.5e-14;			//Gradient energy constant [Jm2/mol]
		//kappa_s1=kappa_s2=kappa_s1/RR/temp/b1/b1;	//Dimensionless
		kappa_s1=kappa_s1/RR/temp/b1/b1;	//Dimensionless
		kappa_s2=kappa_s2/RR/temp/b1/b1;

		//a1_c=b1_c=c1_c=3.6468e-10;					//Lattice constant of Fe (fcc) [nm]
		//atom_n=4.0;  vm0=6.02E23*a1_c*b1_c*c1_c/atom_n;	//Calculation of molar volume (assuming fcc) [m3/mol]

//--- s1 field transformation distortion setting ---
		//eta_s1[1][1]=0.08; eta_s1[2][2]=-0.08; eta_s1[3][3]=0.0;
		//eta_s1[1][1]=-0.06; eta_s1[2][2]=0.06; eta_s1[3][3]=0.0;
		eta_s1[1][1]=eta_s1_11;
		eta_s1[2][2]=eta_s1_22;
		eta_s1[3][3]=eta_s1_33;
		//eta_s1[1][1]=0.04; eta_s1[2][2]=-0.04; eta_s1[3][3]=0.0;
		//eta_s1[1][1]=0.02; eta_s1[2][2]=-0.02; eta_s1[3][3]=0.0;
		//eta_s1[1][2]=eta_s1[2][1]=eta_s1[1][3]=eta_s1[3][1]=eta_s1[2][3]=eta_s1[3][2]=0.0;
		eta_s1[1][2]=eta_s1[2][1]=eta_s1_12;
		eta_s1[1][3]=eta_s1[3][1]=eta_s1_13;
		eta_s1[2][3]=eta_s1[3][2]=eta_s1_23;

//--- s2 field transformation distortion setting ---
		//eta_s2[1][1]=eta_s1[2][2]; 	eta_s2[2][2]=eta_s1[1][1]; 	eta_s2[3][3]=0;
		eta_s2[1][1]=eta_s1[2][2]; 	eta_s2[2][2]=eta_s1[1][1]; 	eta_s2[3][3]=eta_s1[3][3];
		//eta_s2[1][2]=eta_s2[2][1]=eta_s2[1][3]=eta_s2[3][1]=eta_s2[2][3]=eta_s2[3][2]=0.0;
		eta_s2[1][2]=eta_s2[2][1]=eta_s1[1][2];
		eta_s2[1][3]=eta_s2[3][1]=eta_s1[1][3];
		eta_s2[2][3]=eta_s2[3][2]=eta_s1[2][3];

//---Elastic constant (for Fe (fcc))------------------------------------
		//el_fac=1.0E+11*vm0/RR/temp;
		el_fac=el_fac_c*vm0/RR/temp;
		//c11=1.54*el_fac;
		//c44=0.77*el_fac;
		//c12=1.22*el_fac;
		c11=c11c*el_fac;
		c44=c44c*el_fac;
		c12=c12c*el_fac;
		//c12=c11-2.0*c44;
		lam0=c12;	mu0=c44;			//Lame-parameters
		nu0=lam0/2.0/(lam0+mu0);		//Poisson's ratio

//---- External force setting ---------------------------------------------
		//sig22_a=0.0;					//Here, the external force is set to 0.
		//ep11_a=-lam0/4.0/mu0/(lam0+mu0)*sig22_a;	//Assuming plane distortion
		//ep22_a=(lam0+2.0*mu0)/4.0/mu0/(lam0+mu0)*sig22_a;
		//ep12_a=ep21_a=0.0;

//**** Initial field setting of phase field and setting of sin and cos table *******************************************
		prog.ini_field();				//Initial field setting for phase field
		prog.table();					//Setting sin and cos tables and bit transform arrays for Fourier transform

//**** Calculation of time evolution of phase field *******************************************************************************
		while(time1<=time1max){

//---- Display of phase field ----------------------------------------------------
			if((((int)(time1) % 50)==0)){ prog.update_draw(g); }	//Draw a field every multiple of 50 counts
			//if((((int)(time1) % 100)==0)){ prog.repaint(); } 

//---- Save phase field ----------------------------------------------------
			if(time1<=2000.0){Nstep=100;} else{Nstep=500;}
			if((((int)(time1) % Nstep)==0)){ prog.datsave(); }	//Save the field every multiple of Nstep
			//if((((int)(time1) % 200)==0)){ prog.datsave(); }	//Save the field every multiple of 200 counts
			//if(time1==3000.0){ prog.datsave(); }				//Save the field when the count is 3000

//***** Gradient potential *******************************************
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					ip=i+1; im=i-1; jp=j+1; jm=j-1;
					if(i==ndm){ip=0;}  if(i==0){im=ndm;}
					if(j==ndm){jp=0;}  if(j==0){jm=ndm;}
					s1k_su[i][j]=-kappa_s1*(s1h[ip][j]+s1h[im][j]+s1h[i][jp]+s1h[i][jm]
																	-4.0*s1h[i][j]);		//Equation (4)
					s2k_su[i][j]=-kappa_s2*(s2h[ip][j]+s2h[im][j]+s2h[i][jp]+s2h[i][jm]
																	-4.0*s2h[i][j]);		//Equation (4)
					//s1k_su[i][j]=-kappa_s1*( 0.5*(s1h[ip][j]+s1h[im][j]+s1h[i][jp]+s1h[i][jm])
					//												+0.25*(s1h[ip][jp]+s1h[ip][jm]+s1h[im][jp]+s1h[im][jm])
					//												-3.0*s1h[i][j] );		//Equation (4)
					//s2k_su[i][j]=-kappa_s2*( 0.5*(s2h[ip][j]+s2h[im][j]+s2h[i][jp]+s2h[i][jm])
					//												+0.25*(s2h[ip][jp]+s2h[ip][jm]+s2h[im][jp]+s2h[im][jm])
					//												-3.0*s2h[i][j]);		//Equation (4)
				}
			}

//**** Fourier transform of transformation strain field ep11 ********************************
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					xr[i][j]=ep11h0[i][j]=eta_s1[1][1]*s1h[i][j]*s1h[i][j]
															 +eta_s2[1][1]*s2h[i][j]*s2h[i][j];		//Equation (5)
					xi[i][j]=0.0;
				}
			}
			qs=-1.0; prog.rcfft();		//Transform from real space to Fourier space (qs <0)
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){ 
					ep11qrh0[i][j]=xr[i][j];  ep11qih0[i][j]=xi[i][j];
				}
			}
			ep11qrh0[0][0]=ep11qih0[0][0]=0.0;

//**** Fourier transform of transformation strain field ep22 ********************************
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					xr[i][j]=ep22h0[i][j]=eta_s1[2][2]*s1h[i][j]*s1h[i][j]
															 +eta_s2[2][2]*s2h[i][j]*s2h[i][j]; 	//Equation (5)
					xi[i][j]=0.0;
				}
			}
			qs=-1.0; prog.rcfft();		//Transform from real space to Fourier space (qs <0)
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){ 
					ep22qrh0[i][j]=xr[i][j];  ep22qih0[i][j]=xi[i][j];
				}
			}
			ep22qrh0[0][0]=ep22qih0[0][0]=0.0;

//*** Calculation of average value of transformation strain field ***
			sum11=sum22=0.0;
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){ sum11+=ep11h0[i][j];  sum22+=ep22h0[i][j]; }
			}
			ep11_0=sum11/nd/nd;  ep22_0=sum22/nd/nd;

//***** Calculation of constraint strain fluctuation amount ec11 *************************************
			for(i=0;i<=ndm;i++){
				if(i<=nd2-1){ii=i;}  if(i>=nd2){ii=i-nd;}
				for(j=0;j<=ndm;j++){
					if(j<=nd2-1){jj=j;}  if(j>=nd2){jj=j-nd;}
					alnn=Math.sqrt((double)ii*(double)ii+(double)jj*(double)jj);
					if(alnn==0.){alnn=1.;}
					nxx=(double)ii/alnn*(double)ii/alnn;
					nyy=(double)jj/alnn*(double)jj/alnn;
					Z11ep=nxx*(2.0*(1.0-nu0)-nxx-nu0/(1.0-nu0)*nyy)/(1.0-2.0*nu0);
					Z12ep=nxx*(2.0*nu0     -nyy-nu0/(1.0-nu0)*nxx)/(1.0-2.0*nu0);
					xr[i][j]=Z11ep*ep11qrh0[i][j]+Z12ep*ep22qrh0[i][j]; 		//Equation (10)
					xi[i][j]=Z11ep*ep11qih0[i][j]+Z12ep*ep22qih0[i][j]; 		//Equation (10)
				}
			}
			qs=1.0; prog.rcfft();		//Fourier space to real space transformation (qs> 0)
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){ ec11[i][j]=xr[i][j]; }
			}

//***** Calculation of constraint strain fluctuation amount ec22 *****************************
			for(i=0;i<=ndm;i++){
				if(i<=nd2-1){ii=i;}  if(i>=nd2){ii=i-nd;}
				for(j=0;j<=ndm;j++){
					if(j<=nd2-1){jj=j;}  if(j>=nd2){jj=j-nd;}
					alnn=Math.sqrt((double)ii*(double)ii+(double)jj*(double)jj);
					if(alnn==0.){alnn=1.;}
					nxx=(double)ii/alnn*(double)ii/alnn;
					nyy=(double)jj/alnn*(double)jj/alnn;
					Z21ep=nyy*(2.0*nu0     -nxx-nu0/(1.0-nu0)*nyy)/(1.0-2.0*nu0);
					Z22ep=nyy*(2.0*(1.0-nu0)-nyy-nu0/(1.0-nu0)*nxx)/(1.0-2.0*nu0);
					xr[i][j]=Z21ep*ep11qrh0[i][j]+Z22ep*ep22qrh0[i][j]; 		//Equation (10)
					xi[i][j]=Z21ep*ep11qih0[i][j]+Z22ep*ep22qih0[i][j]; 		//Equation (10)
				}
			}
			qs=1.0; prog.rcfft();		//Fourier space to real space transformation (qs> 0)
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){ ec22[i][j]=xr[i][j]; }
			}

//****** Calculation of potential and evolution equations ************************************************************************
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					s1=s1h[i][j];  	s2=s2h[i][j];

//****** Calculation of chemical potential ********************************
					s1k_chem=AA0*(AA1*s1+AA2*s1*s1*s1+AA3*s1*s1*s1*s1*s1
											 +AA4*s1*s2*s2+AA5*s1*(2.0*s1*s1*s2*s2+s2*s2*s2*s2) ); 		//Equation (2)
					s2k_chem=AA0*(AA1*s2+AA2*s2*s2*s2+AA3*s2*s2*s2*s2*s2
											 +AA4*s2*s1*s1+AA5*s2*(2.0*s2*s2*s1*s1+s1*s1*s1*s1) ); 		//Equation (2)

//****** Calculation of elastic potential ********************************
					ep11T=ep11h0[i][j]-ep11_0-ec11[i][j]-ep11_a;
					ep22T=ep22h0[i][j]-ep22_0-ec22[i][j]-ep22_a;

					s1k_str=2.0*s1*( ep11T*((lam0+2.0*mu0)*eta_s1[1][1]+lam0*eta_s1[2][2])
								 					+ep22T*((lam0+2.0*mu0)*eta_s1[2][2]+lam0*eta_s1[1][1]) );	//Equation (7)
					s2k_str=2.0*s2*( ep11T*((lam0+2.0*mu0)*eta_s2[1][1]+lam0*eta_s2[2][2])
								 					+ep22T*((lam0+2.0*mu0)*eta_s2[2][2]+lam0*eta_s2[1][1]) );	//Equation (7)

//****** Calculation of time evolution of phase field ********************************
					//s1ddtt=-smob*(s1k_chem+s1k_su[i][j]); 					//Equation (12)
					//s2ddtt=-smob*(s2k_chem+s2k_su[i][j]); 					//Equation (12)
					s1ddtt=-smob*(s1k_chem+s1k_su[i][j]+s1k_str); 				//Equation (12)
					s2ddtt=-smob*(s2k_chem+s2k_su[i][j]+s2k_str); 				//Equation (12)
					s1h[i][j]=s1h[i][j]+( s1ddtt+ds_fac*(2.0*Math.random()-1.0) )*delt;
					s2h[i][j]=s2h[i][j]+( s2ddtt+ds_fac*(2.0*Math.random()-1.0) )*delt;

//--- Correction of the domain of s (0 <= s <= 1) ---
					if(s1h[i][j]>=1.0){s1h[i][j]=1.0;}  if(s1h[i][j]<=-1.0){s1h[i][j]=-1.0;}
					if(s2h[i][j]>=1.0){s2h[i][j]=1.0;}  if(s2h[i][j]<=-1.0){s2h[i][j]=-1.0;}
				}
			}

//**** [Time increase] *************************************************
			time1=time1+1.0;	//Increased calculation time
		}//while

//----------------------------------------------------------------
		System.out.printf("\n Finished. Click the x in the upper right corner of the graph to finish. \n");
	}//main

// The following is a subroutine.
// **** [Phase field initial settings] *****************************************************
	public void ini_field(){
		int i, j;
		double fac1;

		fac1=0.9; //Maximum initial fluctuation
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				s1h[i][j]=fac1*(2.0*Math.random()-1.0); 
				s2h[i][j]=fac1*(2.0*Math.random()-1.0);		//When placing the nucleus evenly

				//s1h[i][j]=0.0; s2h[i][j]=0.0;
				if(Math.abs(j-nd2)<(nd/40)){
					//s1h[i][j]=2.0*Math.random()-1.0;  s2h[i][j]=2.0*Math.random()-1.0;
				}	//When placing the nucleus in the center
			}
		}

	}

	//******* [Sin and cos tables for Fourier transform and bit transform array settings] **********************
	public void table(){
		int it, it1, it2, mc, mn;
		double q;

		q=2.0*PI/(double)nd;
		for(it=0;it<=nd2-1;it++){ c[it]=Math.cos(q*(double)it); s[it]=Math.sin(q*(double)it); }

		ik[0]=0; mn=nd2; mc=1;
		for(it1=1;it1<=ig;it1++){
			for(it2=0;it2<=mc-1;it2++){ ik[it2+mc]=ik[it2]+mn; }
			mn=mn/2; mc=2*mc;
		}
	}

// *******************************************************************************
	public void update_draw( Graphics g ){ g=getGraphics(); paint(g); }

// **** [Phase field drawing] *******************************************************
	public void paint(Graphics g){
		//g.clearRect(0, 0, width, height);//Clear Window

		int i, j, ii, jj;
		int icol, icol_r, icol_g, icol_b;
		double c1r, c1g, c1b, c2r, c2g, c2b;
		double c_r, c_g, c_b;
		int ixmin=0, iymin=0, igx, igy, irad0;
		int ixmax=xwidth, iymax=yheight;
		double c, x, xmax, xmin, y, ymax, ymin, rad0;

		xmin=0.; xmax=1.; 						//Minimum and maximum values on the horizontal axis
		ymin=0.; ymax=1.; 						//Minimum and maximum values on the vertical axis
		rad0=1.0/(double)nd/2.0; 				//Half the length of the diff block
		irad0=1+(int)( ((double)ixmax-(double)ixmin)/(xmax-xmin)*rad0 ); 	//Pixelization of rad0

		System.out.printf("%f \n", time1); 	//Display the number of calculation repetitions on standard input / output

		for(i=0;i<=nd;i++){
			for(j=0;j<=nd;j++){
				//Position coordinates of phase field (actual value)
				x=1.0/(double)nd*(double)i+rad0;
				y=1.0/(double)nd*(double)j+rad0;
				//Position coordinates of phase field (converted to screen coordinates)
				igx=(int)( ((double)ixmax-(double)ixmin)*(x-xmin)/(xmax-xmin)+(double)ixmin );
				igy=(int)( ((double)iymax-(double)iymin)*(y-ymin)/(ymax-ymin)+(double)iymin );

				//Phase field value of individual diff blocks
				ii=i; jj=j;
				if(i==nd){ii=0;} if(j==nd){jj=0;}			//Periodic boundary conditions

				c1r=c1g=c1b=c2r=c2g=c2b=0.0;
				if(s1h[ii][jj]>=0.0){c1r=s1h[ii][jj];c1g=0.0;c1b=0.0;}
				if(s1h[ii][jj]<0.0) {c1r=0.0;c1g=0.0;c1b=-s1h[ii][jj];}
				if(s2h[ii][jj]>=0.0){c2r=0.0;c2g=s2h[ii][jj];c2b=0.0;}
				if(s2h[ii][jj]<0.0) {c2r=-s2h[ii][jj];c2g=-s2h[ii][jj];c2b=0.0;}
				c_r=c1r+c2r;
				c_g=c1g+c2g;
				c_b=c1b+c2b;
				//c_r= s1h[ii][jj]*s1h[ii][jj];			//s1 red
				//c_g= s2h[ii][jj]*s2h[ii][jj]; 		//s2 green
				//c_b=1.0-c_r-c_g; 						//Blue phase before metamorphosis

				if(c_r>1.0){c_r=1.0;}  if(c_r<0.0){c_r=0.0;}
				if(c_g>1.0){c_g=1.0;}  if(c_g<0.0){c_g=0.0;}
				if(c_b>1.0){c_b=1.0;}  if(c_b<0.0){c_b=0.0;}

				icol_r=(int)(255.*c_r);  	icol_g=(int)(255.*c_g);  icol_b=(int)(255.*c_b);//Convert to 256 layers
				g.setColor(new Color(icol_r, icol_g, icol_b)); 								//Set the color of the difference block
				g.fillRect(insetx+igx-irad0,insety+igy-irad0, irad0*2, irad0*2);			//Individual difference block drawing
			}
		}
	}

//***** [One-dimensional Fast Fourier Transform (FFT)] **************************************
	public void fft(){
		int ix, ka, kb, l2, lf, mf, n2, nf;
		double tj, tr;

		l2=1;
		for(lf=1;lf<=ig;lf++){
			n2=nd2/l2;
			for(mf=1;mf<=l2;mf++){
				for(nf=0;nf<=n2-1;nf++){
					ix=nf*l2; ka=nf+2*n2*(mf-1); kb=ka+n2;
					tr=xrf[ka]-xrf[kb];            tj=xif[ka]-xif[kb];
					xrf[ka]=xrf[ka]+xrf[kb];       xif[ka]=xif[ka]+xif[kb];
					xrf[kb]=tr*c[ix]+tj*qs*s[ix];    xif[kb]=tj*c[ix]-tr*qs*s[ix];
				}
			}
			l2=l2*2;
		}
	}

//**** [Two-dimensional Fast Fourier Transform (RC FFT)] ***********************************
	public void rcfft(){
		int i, ic, ir, j;

		for(ir=0;ir<=ndm;ir++){
			for(ic=0;ic<=ndm;ic++){ xrf[ic]=xr[ir][ic];  xif[ic]=xi[ir][ic]; }
			fft();
			for(ic=0;ic<=ndm;ic++){ xr[ir][ic]=xrf[ik[ic]];  xi[ir][ic]=xif[ik[ic]]; }
		}
		for(ic=0;ic<=ndm;ic++){
			for(ir=0;ir<=ndm;ir++){ xrf[ir]=xr[ir][ic];  xif[ir]=xi[ir][ic]; }
			fft();
			for(ir=0;ir<=ndm;ir++){ xr[ir][ic]=xrf[ik[ir]];  xi[ir][ic]=xif[ik[ir]]; }
		}
		if(qs>0.0){return;}
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){ xr[i][j]=xr[i][j]/nd/nd;  xi[i][j]=xi[i][j]/nd/nd; }
		}
	}

//**** [Data storage] ************************************
	private void datsave() throws Exception{
		int	i, j;

		//Name the save file test.dat.
		PrintWriter outfile= new PrintWriter(				//Open file
			new BufferedWriter(new FileWriter("test.dat", true)) );	//postscript

		outfile.println(time1);								//Write count
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				outfile.println(s1h[i][j]);					//Write phase field
				outfile.println(s2h[i][j]);					//Write phase field
			}
		}
		outfile.close();									//File close
	}

//**** [Read data] ************************************
	private void datin() throws Exception{
		int	i, j;
		String s_data;

		BufferedReader infile=new BufferedReader(new FileReader("ini000.dat"));//Open file

		s_data=infile.readLine();  							//Read as a string
		time1=new Double(s_data).doubleValue();				//Convert letters to numbers
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				s_data=infile.readLine();  					//Read as a string
				s1h[i][j]=new Double(s_data).doubleValue();	//Convert letters to numbers
				s_data=infile.readLine(); 					//Read as a string
				s2h[i][j]=new Double(s_data).doubleValue();	//Convert letters to numbers
			}
		}
		infile.close();										//File close
	}

//*******************************************************************************************************
}//FePt2D_001_v2
//*** The end of the program ********************************************************************************
