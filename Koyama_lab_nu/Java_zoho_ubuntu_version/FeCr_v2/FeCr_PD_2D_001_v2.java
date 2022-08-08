//*** [program (FeCr_PD_2D_001_v2.java)] ************************************************
//*** [import] ****************************
import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.util.*; // for parameters.txt

public class FeCr_PD_2D_001_v2 extends Frame{

//*** [Global variables] *****************************************************************************************
	//static int ND=64;					//Number of divisions on one side of the organization
	//static int nd=ND;					//Number of concentration divisions
	//static int ndm=ND-1;				//Number of concentration divisions-1
	static int ND;						//Number of divisions on one side of the organization
	static int nd;						//Number of concentration divisions
	static int ndm;						//Number of concentration divisions-1
	static int width;					//Overall width of the window
	static int height;					//Overall height of the window
	static int xwidth;					//Drawing area width
	static int yheight;					//Drawing area height
	static int insetx;					//Window frame width (left and right and bottom)
	static int insety;					//Window frame width (top)
	static double PI=3.141592;			//pi
	static double RR=8.3145;			//Gas constant
	//static double [][] ch=new double[ND][ND];	//Concentration data array in tissue
	static double [][] ch;				//Concentration data array in tissue
	static Graphics g;					//Free energy curve screen graphics object
	static double time1;				//Calculation time (count number)
	static double temp; 				//Temperature [K]
	static double c0;					//Alloy composition (mole fraction)


//*** [constructor] ****************************
	public FeCr_PD_2D_001_v2(){
		xwidth=400; yheight=400; 		//Horizontal and vertical length of the drawing screen (in pixels)
		insetx=4; insety=30;			//Frame length of drawing screen
		width=xwidth+insetx*2;  		//Horizontal length of the entire drawing window
		height=yheight+insetx+insety; 	//Vertical length of the entire drawing window
		setSize(width, height);			//Set of drawing windows
		setBackground(Color.white); 	//Set the color of the drawing part of the drawing window to white
		setVisible(true); 				//Make the drawing window visible
		addWindowListener(new WindowAdapter(){ 
			public void windowClosing(WindowEvent e){ System.exit(0); }
										//Operation when closing Window (Setting of x in the upper right of Window)
		});
	}

//*** [main program] *****************************
	public static void main(String[] args) throws Exception{//Exception handling is not performed

		FeCr_PD_2D_001_v2 prog=new FeCr_PD_2D_001_v2();	//Generate an instance prog of FeCr_PD_2D_001_v2

		int i, j; 									//integer
		int ip, im, jp, jm; 						//integer(i+1, i-1, j+1, j-1)
		double delt;								//Time step (dimensionless)
		double al;									//Length of one side of the calculation area
		double time1max;							//Calculation time (maximum count number)
		//double [][] ck = new double[ND][ND];		//Diffusion potential
		//double [][] ch2 = new double[ND][ND]; 		//Concentration data preliminary sequence in tissue
		double [][] ck;								//Diffusion potential
		double [][] ch2;					 		//Concentration data preliminary sequence in tissue
		double mu_chem, mu_surf, mu_str, mu_mag;	//Each potential
		double c1, c2;								//Concentration (Fe: 1, Cr: 2)
		double a01, atomNo;							//Fe (bcc) lattice constant and number of atoms in a unit cell
		double vm0;									//molar volume

		double L0;									//Interatomic interaction parameters
		double L0_0, L0_1;							//Coefficients in interatomic interaction parameters
		double kappa_c;								//Concentration gradient energy coefficient

		double eta;									//Lattice mismatch
		double c11a, c12a, c44a;					//Elastic constant of Fe (bcc)
		double c11b, c12b, c44b;					//Elastic constant of Cr (bcc)
		double c11, c12, c44;						//Overall elastic constant
		double y100;								//Elastic energy function

		//Magnetic excess energy related parameters
		double Tc, d2Tc;							//Curie-temperature and its derivative with respect to its composition
		double tau;									//Temperature standardized by Curie temperature
		double ftau, dftau;							//Magnetic excess energy-function and derivative with respect to temperature
		double Bc, d2Bc;							//Bohr magneton and its derivative with respect to its composition
		double Tca, Tcb, Tcab0, Tcab1;		//Curie temperature-related coefficients
		double Bca, Bcb, Bcab0; 				//Bohr magneton-related coefficients
		double p_mag, D_mag;						//p and D

		double Mc;									//Mobility function and its derivative
		double c_flu;								//The magnitude of the fluctuation of the concentration field
		double cddtt;								//Concentration increment

		double b1;									//Difference block size
		double c2ip, c2im, c2jp, c2jm; 				//In the difference block, centering on c2, the density on the top, bottom, left, and right
		double sumc, dc; 							//Sum of concentration fields, deviation from average composition

		double [] data = new double[40];
		int Nstep;
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
		temp     = data[1];
		c0       = data[2];
		delt     = data[3];
		time1max = data[4];
		al       = data[5];
		vm0      = data[6];
		L0_0     = data[7];
		L0_1     = data[8];
		kappa_c  = data[9];
		eta      = data[10];
		c11a     = data[11];
		c12a     = data[12];
		c44a     = data[13];
		c11b     = data[14];
		c12b     = data[15];
		c44b     = data[16];
		y100     = data[17];
		Tca      = data[18];
		Tcb      = data[19];
		Tcab0    = data[20];
		Tcab1    = data[21];
		Bca      = data[22];
		Bcb      = data[23];
		Bcab0    = data[24];
		p_mag    = data[25];
		D_mag    = data[26];
		Mc       = data[27];
		c_flu    = data[28];
		Nstep    = (int)data[29];
		//
		nd=ND;		//Number of concentration divisions
		ndm=ND-1;	//Number of divisions of concentration-1
		//
		ch  = new double[ND][ND];
		ck  = new double[ND][ND];	//Diffusion potential
		ch2 = new double[ND][ND]; 	//Concentration data preliminary sequence in tissue
		//
		//temp=673.0;  			//Temperature [K]
		//c0=0.5;					//Average composition of alloy (A-50at% B alloy is set here)
		//delt=0.02;				//Time increments
		time1=0.0;				//Number of repetitions of calculation
		//time1max=1.0e+08;		//Maximum number of repetitions

		//al=30.0; 				//Length of one side of 2D calculation area [nm]
		al=al*1.0e-9;			//Convert to [m]
		b1=al/(double)nd;		//Difference 1 block size

		//a01=0.28664e-09;		//Lattice constant
		//atomNo=2.0;				//Number of atoms in a unit cell
		//vm0=6.02E23*a01*a01*a01/atomNo;//Molar volume

		//L0_0=20500.0;			//Interatomic interaction parameter [J/mol] Fe-Cr (bcc)
		//L0_1=-9.68;				//L=L0_0+L0_1*T
		L0=(L0_0+L0_1*temp)/RR/temp; //Dimensionless interatomic interaction parameters

		//kappa_c=2.0e-15;		//Concentration gradient energy coefficient, unit is [Jm2/mol]
		kappa_c=kappa_c/b1/b1/RR/temp;	//Concentration gradient energy coefficient, dimensionless with (b1^2*rr*temp)

		//eta=0.00614;			//Lattice mismatch

		//c11a=2.3310e+11; 		//Elastic modulus of bcc Fe [Pa]
		//c12a=1.3544e+11;
		//c44a=1.1783e+11;

		//c11b=3.500e+11; 		//Elastic constant of bcc Cr [Pa]
		//c12b=0.678e+11;
		//c44b=1.008e+11;

		c11a=c11a*vm0/RR/temp; 		//Dimensionless with RT
		c12a=c12a*vm0/RR/temp;
		c44a=c44a*vm0/RR/temp;
		c11b=c11b*vm0/RR/temp;
		c12b=c12b*vm0/RR/temp;
		c44b=c44b*vm0/RR/temp;

		c11=(1.0-c0)*c11a+c0*c11b;		//Elastic modulus of alloy
		c12=(1.0-c0)*c12a+c0*c12b;
		c44=(1.0-c0)*c44a+c0*c44b;

		//y100=c11+c12-2.0*(c12*c12/c11);	//Modulus function Y <100>
		if(y100==0.0){
			y100=c11+c12-2.0*(c12*c12/c11);	//Modulus function Y <100>
			System.out.printf("Y<100>= %f \n",y100);
		}

		//TcFe=1043.0; TcCr=-311.5; TcbFe0=1650.0; TcbFe1=550.0;//Curie temperature-related coefficients
		//Tca=1043.0;  Tcb=-311.5;  Tcab0=1650.0;  Tcab1=550.0;//Curie temperature-related coefficients
		//BcFe=2.22;   BcCr=-0.008; BcbFe0=-0.85;//Bohr magneton-related coefficients
		//Bca=2.22;    Bcb=-0.008;  Bcab0=-0.85;//Bohr magneton-related coefficients
		//p_mag=0.4;								//(bcc)
		//D_mag=518.0/1125.0+11692.0/15975.0*(1.0/p_mag-1.0);

		//Mc=c0*(1.0-c0);			//Mobility of diffusion
		//c_flu=0.1;

//---- Initial concentration field setting at time 0 ----------------------------------------------------
		prog.ini_comp_field();
		//prog.datin();				//When reading from a file

//---- Calculation of time evolution of concentration field ----------------------------------------------------
		while(time1<=time1max){

//---- Display of concentration field ----------------------------------------------------
			//Concentration field drawing every multiple of 200 counts
			if((((int)(time1) % Nstep)==0)){ prog.update_draw(g); }	//To suppress flicker when drawing
			//if((((int)(time1) % 200)==0)){ prog.repaint(); }

//---- Conservation of concentration field ----------------------------------------------------
			if((((int)(time1) % Nstep)==0)){ prog.datsave(); }	//Concentration field storage every multiple of 200 counts
			//if(time1==3000.0){ prog.datsave(); }				//Concentration field storage when the count number is 3000

//---- Calculation of diffusion potential -----------------------------------------------------
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					ip=i+1; im=i-1; jp=j+1; jm=j-1;
					if(i==ndm){ip=0;}  if(i==0){im=ndm;} 	//Periodic boundary conditions
					if(j==ndm){jp=0;}  if(j==0){jm=ndm;} 	//Periodic boundary conditions

					c2=ch[i][j]; 		c1=1.0-c2; 			//C1 and c2 at position (i, j)
					c2ip=ch[ip][j]; c2im=ch[im][j];
					c2jp=ch[i][jp]; c2jm=ch[i][jm];			//Concentration of c2 before, after, left and right in the difference

		 			mu_chem=L0*(c1-c2)+Math.log(c2)-Math.log(c1); 		//Chemical potential difference
					mu_surf=-2.0*kappa_c*(c2ip+c2im+c2jp+c2jm-4.0*c2);	//Concentration gradient potential
					mu_str=2.0*eta*eta*y100*(c2-c0); 					//Elastic potential

					//Tc=Tca*c1+Tcb*c2;//Curie temperature
					//d2Tc=-Tca+Tcb;//Cr composition derivative of Tc
					//Bc=Bca*c1+Bcb*c2;//Magnetization strength per atom (non-dimensionalized with Bohr magneton)
					//d2Bc=-Bca+Bcb;//Cr composition derivative of Bc

					Tc=Tca*c1+Tcb*c2+c1*c2*(Tcab0+Tcab1*(c2-c1));//Curie temperature
					d2Tc=-Tca+Tcb+2.0*Tcab1*c1*c2+(c1-c2)*(Tcab0+Tcab1*(c2-c1));//Cr composition derivative of Tc
					Bc=Bca*c1+Bcb*c2+Bcab0*c1*c2;//Magnetization strength per atom (non-dimensionalized with Bohr magneton)
					d2Bc=-Bca+Bcb+Bcab0*(c1-c2);//Cr composition derivative of Bc
					if(Tc<0.0){Tc=-Tc;  d2Tc=-d2Tc;}
					if(Bc<0.0){Bc=-Bc;  d2Bc=-d2Bc;}

					tau=temp/Tc; 								//Definition of tau
					if(tau<=1.0){ 
						//Calculation of the tau derivative of f (tau) and f (tau)
						ftau=1.0-1.0/D_mag*(79.0/140.0/p_mag/tau
															+474.0/497.0*(1.0/p_mag-1.0)*( Math.pow(tau,3.0)/6.0
															+Math.pow(tau,9.0)/135.0+Math.pow(tau,15.0)/600.0) );
						dftau=-1.0/D_mag*(-79.0/140.0/p_mag/tau/tau
															+474.0/497.0*(1.0/p_mag-1.0)*( Math.pow(tau,2.0)/2.0
															+Math.pow(tau,8.0)/15.0+Math.pow(tau,14.0)/40.0) );
					}
					else{	
						ftau=-1.0/D_mag*(Math.pow(tau,-5.0)/10.0+Math.pow(tau,-15.0)/315.0
														+Math.pow(tau,-25.0)/1500.0);
						dftau=1.0/D_mag*(Math.pow(tau,-6.0)/2.0+Math.pow(tau,-16.0)/21.0
														+Math.pow(tau,-26.0)/60.0);
					}

					mu_mag=ftau/(Bc+1.0)*d2Bc-dftau*tau/Tc*d2Tc*Math.log(Bc+1.0); 
																//Magnetic excess energy potential

					ck[i][j]=mu_chem+mu_surf+mu_str+mu_mag;		//Diffusion potential

				}
			}

//---- Time change of concentration field (differential decomposition explicit method of nonlinear diffusion equation) ------------------------------------
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					ip=i+1; im=i-1;  jp=j+1;  jm=j-1;
					if(i==ndm) {ip=0;}	if(i==0) {im=ndm;} 	//Periodic boundary conditions
					if(j==ndm) {jp=0;}	if(j==0) {jm=ndm;} 	//Periodic boundary conditions
					cddtt=Mc*(ck[ip][j]+ck[im][j]+ck[i][jp]+ck[i][jm]-4.0*ck[i][j]);	//Non-linear diffusion equation
					//ch2[i][j]=ch[i][j]+cddtt*delt; 									//Time evolution of concentration field
					ch2[i][j]=ch[i][j]+( cddtt+c_flu*(2.0*Math.random()-1.0) ) *delt;	//Time evolution of concentration field (introduction of concentration fluctuation)
				}
			}

//*** [Correction of the balance of the concentration field] *******************************************************
//*** Since it is a numerical calculation, the balance of the concentration field is corrected (actually, it is not necessary to perform each step).] ****
  		sumc=0.;
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					sumc=sumc+ch2[i][j]; 					//Integral of concentration field
				}
			}
			dc=sumc/(double)nd/(double)nd-c0;				//Fluctuation amount of concentration field

			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					ch[i][j]=ch2[i][j]-dc; 					//Concentration field correction
					if(ch[i][j]>=1.){ch[i][j]=1.0-1.0e-6;} 	//Correction when the density exceeds 1
					if(ch[i][j]<=0.){ch[i][j]=1.0e-6;} 		//Correction when the density is less than 0
				}
			}

//******[Time increase]*************************************************
			time1=time1+1.0; 	//Increased calculation time
		}//while

//----------------------------------------------------------------

		System.out.printf("\n Finished. Click the x in the upper right corner of the graph to finish. \n");

}//main

// The following is a subroutine.
// **** Initial concentration field setting *****************************************************
	public void ini_comp_field(){
		int i, j;
		double rnd0, fac1;

		fac1=0.01; 	//Set the maximum amount of change in initial concentration fluctuation to 1%
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				rnd0=2.0*Math.random()-1.0;  ch[i][j]=c0+rnd0*fac1; 	//Set the initial concentration field with random numbers
			}
		}

	}

// *********************************************************
//Separately define draw function to prevent flicker
	public void update_draw( Graphics g ){ g=getGraphics(); paint(g); }

// **** Density field drawing *****************************************************
	public void paint(Graphics g){
		//g.clearRect(0, 0, width, height);		//Clear Window

		int i, j, ii, jj;
		double c, x, xmax, xmin, y, ymax, ymin, rad0;
		int ixmin=0, iymin=0, igx, igy, irad0;
		int ixmax=xwidth, iymax=yheight;
		int icol;

		xmin=0.; xmax=1.; 						//Minimum and maximum values on the horizontal axis
		ymin=0.; ymax=1.; 						//Minimum and maximum values on the vertical axis

		rad0=1.0/(double)nd/2.0; 				//Half the length of the diff block
		irad0=1+(int)( ((double)ixmax-(double)ixmin)/(xmax-xmin)*rad0 ); 	//Pixelization of rad0

		System.out.printf("%f \n", time1); 	//Display the number of calculation repetitions on standard input / output

		for(i=0;i<=nd;i++){
			for(j=0;j<=nd;j++){
				//Position coordinates of concentration field (actual value)
				x=1.0/(double)nd*(double)i+rad0;
				y=1.0/(double)nd*(double)j+rad0;
				//Position coordinates of concentration field (converted to screen coordinates)
				igx=(int)( ((double)ixmax-(double)ixmin)*(x-xmin)/(xmax-xmin)+(double)ixmin );
				igy=(int)( ((double)iymax-(double)iymin)*(y-ymin)/(ymax-ymin)+(double)iymin );

				//Concentration value of individual difference blocks
				ii=i; jj=j;
				if(i==nd){ii=0;} if(j==nd){jj=0;}								//Periodic boundary conditions
				icol=(int)(255.0*(1.0-ch[ii][jj]));								//Make the color gradation grayscale
				//icol=(int)(255.0*ch[ii][jj]);									//When reversing light and dark
				if(icol>=255){icol=255;} if(icol<=0){icol=0;}					//Correction of light and dark range
				g.setColor(new Color(icol,icol,icol)); 							//Set the density in light and dark
				g.fillRect(insetx+igx-irad0,insety+igy-irad0, irad0*2, irad0*2);//Individual difference block drawing
			}
		}
	}

//*** [Data storage] ************************************
	private void datsave() throws Exception{
		int	i, j;

		//Name the save file test.dat.
		PrintWriter outfile= new PrintWriter(
			new BufferedWriter(new FileWriter("test.dat", true)) );	//File open postscript

		outfile.println(time1);					//Write count
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				outfile.println(ch[i][j]);		//Writing the concentration field
			}
		}
		outfile.close();						//File close
	}

//*** [Read data] ************************************
	private void datin() throws Exception{
		int	i, j;
		String s_data;

		BufferedReader infile=new BufferedReader(new FileReader("ini000.dat"));//Open file

		s_data=infile.readLine();  									//Read as a character string
		time1=new Double(s_data).doubleValue();						//Convert letters to numbers
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				s_data=infile.readLine();  							//Read as a character string
				ch[i][j]=new Double(s_data).doubleValue();			//Convert letters to numbers
			}
		}
		infile.close();												//File close
	}

//****************************************************************
}//FeCr_PD_2D_001_v2
//*** The end of the program ************************************************************
