//*** [program (PD2D_001_v2.java)] ************************************************
//*** [import] ****************************
import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.util.*; // for parameters.txt

public class PD2D_001_v2 extends Frame{

//*** [Global variables] *****************************************************************************************
	static int ND=64;			//Number of divisions on one side of the organization
	static int nd=ND;			//Number of concentration divisions
	static int ndm=ND-1;		//Number of divisions of concentration-1
	static int width;			//Overall width of the window
	static int height;			//Overall height of the window
	static int xwidth;			//Drawing area width
	static int yheight;			//Drawing area height
	static int insetx;			//Window frame width (left and right and bottom)
	static int insety;			//Window frame width (top)
	static double PI=3.141592;	//pi
	static double RR=8.3145;	//Gas constant
	static double [][] ch=new double[ND][ND];	//Concentration data array in tissue
	static Graphics g;			//Free energy curve screen graphics object
	static double time1;		//Calculation time (count number)
	static double temp; 		//Temperature [K]
	static double c0;			//Alloy composition (mole fraction)


//*** [constructor] ****************************
	public PD2D_001_v2(){
		xwidth=400; yheight=400;		//Horizontal and vertical length of the drawing screen (in pixels)
		insetx=4; insety=30;			//Frame length of drawing screen
		width=xwidth+insetx*2;  		//Horizontal length of the entire drawing window
		height=yheight+insetx+insety;	//Vertical length of the entire drawing window
		setSize(width, height);			//Set of drawing windows
		setBackground(Color.white); 	//Set the color of the drawing part of the drawing window to white
		setVisible(true); 				//Make the drawing window visible
		addWindowListener(new WindowAdapter(){ 
			public void windowClosing(WindowEvent e){ System.exit(0); }
										//Operation when closing indow (setting of x in the upper right of Window)
		});
	}

//*** [main program] *****************************
	public static void main(String[] args) throws Exception{//Exception handling is not performed

		PD2D_001_v2 prog=new PD2D_001_v2();			//Generate instance prog of PD2D_001_v2

		int i, j; 								//integer
		int ip, im, jp, jm; 					//integer (i+1, i-1, j+1, j-1)
		double delt;							//Time step (dimensionless)
		double al;								//Length of one side of the calculation area
		double time1max;						//Calculation time (maximum count number)
		double [][] ck = new double[ND][ND];	//Diffusion potential
		double [][] ch2 = new double[ND][ND]; 	//Concentration data preliminary sequence in tissue
		double mu_chem, mu_surf;				//Each potential
		double c1, c2;							//Concentration (Fe: 1, Cu: 2)
		double L0;								//Interatomic interaction parameters
		double kappa_c;							//Concentration gradient energy coefficient
		double Mc, M0;							//Mobility function and its derivative
		double c_flu;							//The magnitude of the fluctuation of the concentration field
		double cddtt;							//Concentration increment

		double b1;								//Difference block size
		double c2ip, c2im, c2jp, c2jm; 			//In the difference block, centering on c2, the density on the top, bottom, left, and right
		double sumc, dc; 						//Sum of concentration fields, deviation from average composition
		
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
		temp     = data[1];
		c0       = data[2];
		delt     = data[3];
		time1max = data[4];
		al       = data[5];
		L0       = data[6];
		kappa_c  = data[7];
		M0       = data[8];
		c_flu    = data[9];
		//
		//temp=1000.0;  			//[K]
		//c0=0.4;					//Average composition of alloy (A-40at% B alloy is set here)
		//delt=0.04;				//Time increments
		time1=0.0;				//Number of repetitions of calculation
		//time1max=1.0e+08;		//Maximum number of repetitions

		//al=60.0; 				//Length of one side of 2D calculation area [nm]
		al=al*1.0e-9;			//Convert to [m]
		b1=al/(double)nd;		//Difference 1 block size

		//L0=2.5e+04;				//Interatomic interaction parameter [J/mol]
		L0=L0/RR/temp;			//Dimensionless

		//kappa_c=5.0e-15;		//Concentration gradient energy coefficient, unit is [Jm2/mol]
		kappa_c=kappa_c/b1/b1/RR/temp;	//Concentration gradient energy coefficient, dimensionless with (b1^2*rr*temp)

		//Mc=c0*(1.0-c0);			//Mobility of diffusion
		Mc=c0*(M0-c0);			//Mobility of diffusion
		c_flu=0.1;

//---- Initial concentration field setting at time 0 ----------------------------------------------------
		prog.ini_comp_field();
		//prog.datin();				//When reading from a file

//---- Calculation of time evolution of concentration field ----------------------------------------------------
		while(time1<=time1max){

//---- Display of concentration field ----------------------------------------------------
			//Concentration field drawing every multiple of 200 counts
			if((((int)(time1) % 200)==0)){ prog.update_draw(g); }	//To suppress flicker when drawing
			//if((((int)(time1) % 200)==0)){ prog.repaint(); }

//---- Conservation of concentration field ----------------------------------------------------
			if((((int)(time1) % 500)==0)){ prog.datsave(); }		//Concentration field storage every multiple of 500 counts
			//if(time1==3000.0){ prog.datsave(); }					//Concentration field storage when the count number is 3000

//---- Calculation of diffusion potential -----------------------------------------------------
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					ip=i+1; im=i-1; jp=j+1; jm=j-1;
					if(i==ndm){ip=0;}  if(i==0){im=ndm;} 	//Periodic boundary conditions
					if(j==ndm){jp=0;}  if(j==0){jm=ndm;} 	//Periodic boundary conditions

					c2=ch[i][j]; 		c1=1.0-c2; 			//C1 and c2 at position (i, j)
					c2ip=ch[ip][j]; c2im=ch[im][j]; c2jp=ch[i][jp]; c2jm=ch[i][jm]; //Concentration of c2 before, after, left and right in the difference

		 			mu_chem=L0*(c1-c2)+Math.log(c2)-Math.log(c1); 		//Chemical potential difference
					mu_surf=-2.0*kappa_c*(c2ip+c2im+c2jp+c2jm-4.0*c2);	//Concentration gradient potential
					ck[i][j]=mu_chem+mu_surf; 							//Diffusion potential
				}
			}

//---- Time change of concentration field (differential decomposition explicit method of nonlinear diffusion equation) ------------------------------------
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					ip=i+1; im=i-1;  jp=j+1;  jm=j-1;
					if(i==ndm) {ip=0;}	if(i==0) {im=ndm;} 	//Periodic boundary conditions
					if(j==ndm) {jp=0;}	if(j==0) {jm=ndm;} 	//Periodic boundary conditions
					cddtt=Mc*(ck[ip][j]+ck[im][j]+ck[i][jp]+ck[i][jm]-4.0* ck[i][j]);	//Non-linear diffusion equation
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
				rnd0=2.0*Math.random()-1.0;  ch[i][j]=c0+rnd0*fac1;	//Set the initial concentration field with random numbers
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
		irad0=1+(int)( ((double)ixmax-(double)ixmin)/(xmax-xmin)*rad0 );	//Pixelization of rad0

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
				if(i==nd){ii=0;} if(j==nd){jj=0;}									//Periodic boundary conditions
				icol=(int)(255.0*(1.0-ch[ii][jj]));									//Make the color gradation grayscale
				//icol=(int)(255.0*ch[ii][jj]);										//When reversing light and dark
				if(icol>=255){icol=255;} if(icol<=0){icol=0;}						//Correction of light and dark range
				g.setColor(new Color(icol,icol,icol)); 								//Set the density in light and dark
				g.fillRect(insetx+igx-irad0,insety+igy-irad0, irad0*2, irad0*2);	//Individual difference block drawing
			}
		}
	}

//*** [Data storage] ************************************
	private void datsave() throws Exception{
		int	i, j;

		//Name the save file test.dat
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

		s_data=infile.readLine();  							//Read as a character string
		time1=new Double(s_data).doubleValue();				//Convert letters to numbers
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				s_data=infile.readLine();  					//Read as a character string
				ch[i][j]=new Double(s_data).doubleValue();	//Convert letters to numbers
			}
		}
		infile.close();										//File close
	}

//****************************************************************
}//PD2D_001_v2
//*** The end of the program ************************************************************
