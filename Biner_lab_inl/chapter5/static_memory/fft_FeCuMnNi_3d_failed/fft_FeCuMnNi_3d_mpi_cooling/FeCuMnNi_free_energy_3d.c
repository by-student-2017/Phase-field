#include <math.h> //M_PI
#include <stdlib.h> //rand()

/* Note
  Although we are using an array of complex numbers,
   the imaginary component of the input and output in this function is 0.
  
  This function needs to be rewritten depending on the system.
  
  If the systems are different, you will rewrite the free energy and
   its derivative obtained from CALPHAD.
  
  I have included comments so that if you understand it while using it in 
   conjunction with the textbook, it will not be difficult to rewrite it.
  
  Good luck !
*/

/* Variable and array list
  Nx: Number of grid points in the x-direction
  Ny: Number of grid points in the y-direction
  Nz: Number of grid points in the z-direction
  tempr: Temperature [K]
  cu(Nx,Ny,Nz): Modulated Cu values at the grid points
  mn(Nx,Ny,Nz): Modulated Mn values at the grid points
  ni(Nx,Ny,Nz): Modulated Ni values at the grid points
  orp(Nx,Ny,Nz): Modulated order parameter values at the grid points
  dgdcu(Nx,Ny,Nz): Functional derivative with respect to Cu
  dgdmn(Nx,Ny,Nz): Functional derivative with respect to Mn
  dgdni(Nx,Ny,Nz): Functional derivative with respect to Ni
  dgdor(Nx,Ny,Nz): Functional derivative with respect to
    non-conserved order parameter */

void FeCuMnNi_free_energy_3d(int Nx, int Ny, int Nz, 
	double *cu, double *mn, double *ni, double *orp, double tempr,
	double *dgdcu, double *dgdmn, double *dgdni, double *dgdor){
	
	double R=8.314472; //The value of gas constant, R [J/mol/K]
	double T=tempr; //Temperature [K]
	double RT=R*T; //RT [J/mol], 6842.810456 [J/mol] for 823.0 [K]
	double beta=1.0/RT; //inverse temperature, 1.4613878411949395e-4 for 823.0 [K]
	
	double constw=5.0e3/RT; //The value W in Eq.5.25(Total free energy, F). The W is [J/mol] unit
		/* The term W*g(eta) in the integrand represents the energy barrier for
		   the phase transformation between the alpha and gamma phases. g(eta)=eta*(1-eta) */
	double Vm=7.09e-6; //Molar volume [m^3/mol]
	double Es=214.0e9; //Energy function of elastic stifness, nGPa
	double conste=Es*Vm/RT; //The normalized value of Y in Eq.5.25(Total free energy, F)
	
	//Lattice mismatch (=lattice misfit)
	//Misfit strain values in Eq.5.26 for each alloying element
	double eta2=3.29e-2;
	double eta3=5.22e-4;
	double eta4=4.75e-4;
	
	//Initial concentrations of alloying elements
	double c02=0.15; //average Cu concentration
	double c03=0.01; //average Mn concentration
	double c04=0.01; //average Ni concentration
	
	//Value of coefficient which will be used as a penalty parameter later
	double Coef=10.0;
	
	//----- ----- ----- ----- ----- ----- ----- -----
	// alpha phase
	//Gibbs energy (Memo: ln(T)=log(T) on C language)
	double G0A1 = 0.0;
	double G0A2 =   4017.0 + (-1.255*T); 
		//2984.135 for 823.0 [K]
	double G0A3 =  -3235.3 +  127.85*T + (-23.7*T*log(T)) + (-0.00742717*T*T) + 60000.0/T;
		//-33909.365609 for 823.0 [K]
	double G0A4 = 8175.084 + (-3.556*T);
		//5248.496 for 823.0 [K]
	//
	double LA12 =  41033.0 + (-6.022*T); 
		//36076.894 for 823.0 [K]
	double LA13 =  -2759.0 +   1.237*T;
		//-1740.949 for 823.0 [K]
	double LA14p=(1789.03-1.92912*T); 
	double LA14;
	//double LA14 = -956.63-1.28726*T + LA14p*(c1-c4);
		//-2016.04498 + 201.36424*(c1-c4) for 823.0 [K]
	double LA23p=-9865.0;
	double LA23;
	//double LA23 =     11190.0-6.0*T + LA23p*(c2-c3);
		//6252.0 - 9865.0*(c2-c3) for 823.0 [K]
	double LA24 =     8366.0+2.80*T;
		//10670.4 for 823.0 [K]
	double LA34p=6276.0;
	double LA34;
	//double LA34 =  -51638.31+3.64*T + LA34p*(c3-c4);
		//-48642.59 + 6276.0*(c3-c4) for 823.0 [K]
	//
	double LA123=30000.0;
	//double LA124=0.0;
	//double LA134=0.0;
	//double LA234=0.0;
	//----- ----- ----- ----- ----- ----- ----- -----
	// gamma phase
	//Gibbs energy (Memo: ln(T)=log(T) on C language)
	double G0G1 =      -1462.4+8.282*T +    (-1.15*T*log(T)) +  6.4e-4*T*T;
		//3027.89638596 for 823.0 [K]
	double G0G2 = 0.0;
	double G0G3 =    -3439.3+131.884*T + (-24.5177*T*log(T)) + (-0.006*T*T) + 69600.0/T;
		//42294.693153 got 823.0 [K]
	double G0G4 = 0.0;
	//
	double LG12p=(11512.0-7.095*T);
	double LG12;
	//double LG12 =     53360.0-12.626*T + LG12p*(c2-c1);
		//42968.802 + 5672.815*(c2-c1) for 823.0 [K]
	double LG13p=-259.0;
	double LG13;
	//double LG13 =        -7762+3.865*T + LG13p*(c1-c3);
		//-4581.105 - 259.0*(c1-c3) for 823.0 [K]
	double LG14p1=(11082.1315-4.45077*T);
	double LG14p2=-725.8051;
	double LG14;
	//double LG14 = -12054.335+3.27413*T + LG14p1*(c1-c4) + LG14p2*(c1-c4)*(c1-c4);
		//-9359.72601 + 7419.14779*(c1-c4) -725.8051*(c1-c4)*(c1-c4) for 823.0 [K]
	double LG23p1=(-10600.0+3.0*T);
	double LG23p2=(-4850.0+3.5*T);
	double LG23;
	//double LG23 =        11820.0-2.3*T + LG23p1*(c2-c3) + LG23p2*(c2-c3)*(c2-c3)*(c2-c3);
		//9927.1 + -8131.0*(c2-c3) + -1969.5*(c2-c3)*(c2-c3)*(c2-c3) for 823.0 [K]
	double LG24p=(-4359.6+1.812*T);
	double LG24;
	//double LG24 =       8366.0+2.802*T + LG24p*(c2-c4);
		//10672.046 + -2868.324*(c2-c4) for 823.0 [K]
	double LG34p=6276.0;
	double LG34;
	//double LG34 =    -58158.0+10.878*T + LG34p*(c3-c4);
		//-49205.406 + 6276*(c3-c4) for 823.0 [K]
	//
	double LG123 =     -68000.0+50.0*T;
		//-26850.0 for 823.0 [K]
	double LG124 =     -73272.0+30.9*T;
		//-47841.3 for 823.0 [K]
	//----- ----- ----- ----- ----- ----- ----- -----	
	int ii;
	
	double c1,c2,c3,c4;
	double funch; //h(eta)
	double funcg; //g(eta)
	double cc0eta; //sum(etai*(ci-c0i))
	double elaste;
	double delascu, delasmn, delasni;
	double dgcua, dgmna, dgnia; // dG/dc for alpha phase
	double dgcug, dgmng, dgnig; // dG/dc for gamma phase
	double gcuni;
	double gcunia; // sum(Ga(ci)) of alpha phase
	double gcunig; // sum(Gg(ci)) of gamma phase
	double thrc=0.0;
	double dgcu;
	//double dhdor;
	
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			for(int k=0;k<Nz;k++){
				ii=i*Ny*Nz+j*Nz+k;
				
				//Current values of the alloying elements at the curent grid point
				c2=cu[ii]; // concentration of Cu(x,y,z)
				c3=mn[ii]; // concentration of Mn(x,y,z)
				c4=ni[ii]; // concentration of Ni(x,y,z)
				c1=1.0-c2-c3-c4; // concentration of Fe(x,y,z)
				
				//Evaluate the value of function h(eta)=eta*eta*(3-2*eta), eta=orp[ii][0]
				funch=orp[ii]*orp[ii]*(3.0-2.0*orp[ii]);
				
				//Evaluate the value of function g(eta)=eta*(1-eta), eta=orp[ii][0]
				funcg=orp[ii]*(1.0-orp[ii]);
				
				//This is related with Eq.5.26
				cc0eta=( eta2*(c2-c02) + eta3*(c3-c03) + eta4*(c4-c04) );
				
				//Calculate the elastic energy
				elaste=conste*cc0eta*cc0eta; //Y*Vm*epsilon0(ci)^2 in Eq.5.25
				
				//Derivative of elastic energy with respect to each alloying element
				/* The quantitiy Y*Vm*epsilon0(ci)^2 is the elastic strain energy
				   include from the coherent phase separation, 
				   where Y is the average stifness and Vm is the molar volume and
				   epsilon0(ci) is the eigenstrain and taken as Eq.5.26 */
				delascu=conste*2.0*cc0eta*eta2; // d(elaste)/d(c2)
				delasmn=conste*2.0*cc0eta*eta3; // d(elaste)/d(c3)
				delasni=conste*2.0*cc0eta*eta4; // d(elaste)/d(c4)
				
				//parameters
				//----- ----- ----- -----
				// alpha phase
				//Gibbs energy (Memo: ln(T)=log(T) on C language)
				G0A1 =      0.0;
				G0A2 =   4017.0 + (-1.255*T); 
				G0A3 =  -3235.3 +  127.85*T + (-23.7*T*log(T)) + (-0.00742717*T*T) + 60000.0/T;
				G0A4 = 8175.084 + (-3.556*T);
				//
				LA12 =  41033.0 + (-6.022*T); 
				LA13 =  -2759.0 +   1.237*T;
				LA14p=(1789.03-1.92912*T); 
				LA14 =    -956.63-1.28726*T + LA14p*(c1-c4);
				LA23p=-9865.0;
				LA23 =        11190.0-6.0*T + LA23p*(c2-c3);
				LA24 =        8366.0+2.80*T;
				LA34p=6276.0;
				LA34 =     -51638.31+3.64*T + LA34p*(c3-c4);
				//
				LA123=  30000.0;
				//LA124=      0.0;
				//LA134=      0.0;
				//LA234=      0.0;
				//----- ----- ----- -----
				// gamma phase
				//Gibbs energy (Memo: ln(T)=log(T) on C language)
				G0G1 =      -1462.4+8.282*T +    (-1.15*T*log(T)) +  6.4e-4*T*T;
				G0G2 =      0.0;
				G0G3 =    -3439.3+131.884*T + (-24.5177*T*log(T)) + (-0.006*T*T) + 69600.0/T;
				G0G4 =      0.0;
				//
				LG12p=(11512.0-7.095*T);
				LG12 =     53360.0-12.626*T + LG12p*(c2-c1);
				LG13p=-259.0;
				LG13 =        -7762+3.865*T + LG13p*(c1-c3);
				LG14p1=(11082.1315-4.45077*T);
				LG14p2=-725.8051;
				LG14 = -12054.335+3.27413*T + LG14p1*(c1-c4) + LG14p2*(c1-c4)*(c1-c4);
				LG23p1=(-10600.0+3.0*T);
				LG23p2=(-4850.0+3.5*T);
				LG23 =        11820.0-2.3*T + LG23p1*(c2-c3) + LG23p2*(c2-c3)*(c2-c3)*(c2-c3);
				LG24 =       8366.0+2.802*T + LG24p*(c2-c4);
				LG34p=6276.0;
				LG34 =    -58158.0+10.878*T + LG34p*(c3-c4);
				//
				LG123 =     -68000.0+50.0*T;
				LG124 =     -73272.0+30.9*T;
				//----- ----- ----- -----
				
				//If Fe concentration is within the limits continue
				if(c1<0.9995 && c1>0.0005){
					//
					// 1:Fe, 2:Cu, 3:Mn, 4:Ni
					
					// Calculate the functional derivatives for alpha phase
					//
					// alpha phase
					//----- ----- ----- ----- ----- -----
					//Functional derivative with respect to Cu concentration
					//----- -----
					//dGa/dCcu=dGa/dc2 (modified original code)
					//d(LA14)*c1*c4=-LA14p*c1*c4
					//d(LA23)*c2*c3=LA23p*c2*c3
					//=G0A2-G0A1
					//+(-LA14p*c1*c4-LA12*c2-LA13*c3-LA14*c4 +LA12*c1 +LA23p*c2*c3+LA23*c3 +LA24*c4)
					//-LA123*c2*c3+LA123*c1*c3
					//+Gmg(omitted)+-RT*log(c1)-RT+RT*log(c2)+RT
					//Gmg:Ref https://doi.org/10.2320/matertrans.47.2765
					dgcua=beta*
						(G0A2-G0A1
						-LA14p*c1*c4
						-LA12*c2 -LA13*c3 -LA14*c4
						+LA12*c1
						+LA23p*c2*c3
						+LA23*c3
						+LA24*c4
						-LA123*c2*c3
						+LA123*c1*c3
						+RT*(log(c2)-log(c1))
						);
					//----- ----- ----- ----- ----- -----
					//Functional derivative with respect to Mn concentration
					//dGa/dCmn=dGa/dc3 (modified original code)
					//d(LA14)*c1*c4=-LA14p*c1*c4
					//d(LA23)*c2*c3=-LA23p*c2*c3
					//d(LA34)*c3*c4=LA34p*c3*c4
					//=G0A3-G0A1
					//+(-LA14p*c1*c4-LA12*c2-LA13*c3-LA14*c4 +LA13*c1 -LA23p*c2*c3+LA23*c2 +LA34p*c3*c4+LA34*c4)
					//-LA123*c2*c3+LA123*c1*c2
					//+Gmg(omitted)+-RT*log(c1)-RT+RT*log(c3)+RT
						// LA123*c1*c2*c3=LA123*c1*c2-LA123*c2*c3
					dgmna=beta*
						(G0A3-G0A1
						-LA14p*c1*c4
						-LA12*c2 -LA13*c3 -LA14*c4
						+LA13*c1
						-LA23p*c2*c3
						+LA23*c2
						+LA34p*c3*c4 +LA34*c4
						-LA123*c2*c3
						+LA123*c1*c2
						+RT*(log(c3)-log(c1))
						);
					//----- ----- ----- ----- ----- -----
					//Functional derivative with respect to Ni concentration
					//dGa/dCni=dGa/dc4 (modified original code)
					//d(LA14)*c1*c4=LA14p*-2*c1*c4
					//d(LA34)*c3*c4=-LA34p*c3*c4
					//=G0A4-G0A1
					//+(LA14p*-2*c1*c4-LA12*c2-LA13*c3-LA14*c4 +LA14*c1 +LA24*c2 -LA34p*c3*c4+LA34*c3)
					//-LA123*c2*c3
					//+Gmg(omitted)+-RT*log(c1)-RT+RT*log(c4)+RT
					dgnia=beta*
						(G0A4-G0A1
						+LA14p*-2.0*c1*c4
						-LA12*c2
						-LA13*c3
						-LA14*c4
						+LA14*c1
						+LA24*c2
						-LA34p*c3*c4 +LA34*c3
						-LA123*c2*c3
						+RT*(log(c4)-log(c1))
						);
					//----- ----- ----- ----- ----- -----
					
					// Calculate the functional derivatives for gamma phase
					//
					// gamma phase
					//----- ----- ----- ----- ----- -----
					//Functional derivative with respect to Cu concentration
					//dGg/dCcu=dGg/dc2 (modified original code)
					//d(LG12)*c1*c2=LG12p*2.0*c1*c2
					//d(LG13)*c1*c3=-LG13p*c1*c3
					//d(LG14)*c1*c4=-LG14p1*c1*c4-LG14p2*2.0*(c1-c4)*c1*c4
					//d(LG23)*c2*c3=LG23p1*c2*c3+LG23p2*3.0*(c2-c3)*(c2-c3)*c2*c3
					//d(LG24)*c2*c4=LG24p*c2*c4
					//=G0G2-G0G1
					//+(LG12p*2.0*c1*c2-LG13p*c1*c3-LG14p1*c1*c4-LG14p2*2.0*(c1-c4)*c1*c4
					//-LG12*c2-LG13*c3-LG14*c4 +LG12*c1 +LG23*c3 +LG24*c4
					//+LG23p1*c2*c3+LG23p2*3.0*(c2-c3)*(c2-c3)*c2*c3+LG24p*c2*c4)
					//-LG123*c2*c3+LG123*c1*c3 -LG124*c2*c4+LG124*c1*c4
					//+Gmg(omitted)+-RT*log(c1)-RT+RT*log(c2)+RT
					//Gmg:Ref https://doi.org/10.2320/matertrans.47.2765
					dgcug=beta*
						(G0G2-G0G1 //+566.3008361308123 ?
						+LG12p*2.0*c1*c2
						-LG13p*c1*c3 //a part of -26591.0*c1*c3 = (-LG13p+LG123)*c1*c3
						-LG14p1*c1*c4-LG14p2*2.0*(c1-c4)*c1*c4
						-LG12*c2
						-LG13*c3
						-LG14*c4
						+LG12*c1
						+LG23p1*c2*c3+LG23p2*3.0*(c2-c3)*(c2-c3)*c2*c3
						+LG23*c3
						+LG24p*c2*c4
						+LG24*c4
						-LG124*c2*c4
						+LG124*c1*c4
						-LG123*c2*c3
						+LG123*c1*c3
						+RT*(log(c2)-log(c1))
						);
					//----- ----- ----- ----- ----- -----
					//Functional derivative with respect to Mn concentration
					//dGg/dCmn=dGg/dc3 (modified original code)
					//d(LG12)*c1*c2=LG12p*c1*c2
					//d(LG13)*c1*c3=LG13p*-2.0*c1*c3
					//d(LG14)*c1*c4=-LG14p1*c1*c4 +LG14p2*-2.0*(c1-c4)*c1*c4
					//d(LG23)*c2*c3=LG23p1*c2*c3 +LG23p2*3.0*(c2-c3)*(c2-c3)*c2*c3
					//d(LG34)*c3*c4=LG34p*c3*c4
					//=G0G3-G0G1
					//+LG12p*c1*c2
					//+LG13p*-2.0*c1*c3
					//-LG14p1*c1*c4 +LG14p2*-2.0*(c1-c4)*c1*c4
					//-LG12*c2-LG13*c3-LG14*c4
					//-LG23p1*c2*c3 +LG23p2*-3.0*(c2-c3)*(c2-c3)*c2*c3 +LG23*c2
					//+LG13*c1
					//+LG34p*c3*c4 +LG34*c4
					//-LG123*c2*c3+LG123*c1*c2 -LG124*c2*c4
					//+Gmg(omitted)+-RT*log(c1)-RT+RT*log(c3)+RT
					dgmng=beta*
						(G0G3-G0G1
						+LG12p*c1*c2
						+LG13p*-2.0*c1*c3
						-LG14p1*c1*c4 +LG14p2*-2.0*(c1-c4)*c1*c4
						-LG12*c2
						-LG13*c3
						-LG14*c4
						-LG23p1*c2*c3 +LG23p2*-3.0*(c2-c3)*(c2-c3)*c2*c3
						+LG23*c2
						+LG13*c1
						+LG34p*c3*c4
						+LG34*c4
						-LG123*c2*c3
						+LG123*c1*c2
						-LG124*c2*c4
						+RT*(log(c3)-log(c1))
						);
					//----- ----- ----- ----- ----- -----
					//Functional derivative with respect to Ni concentration
					//dGg/dCni=dGg/dc4 (modified original code)
					//d(LG12)*c1*c2=LG12p*c1*c2
					//d(LG13)*c1*c3=-LG13p*c1*c3
					//d(LG14)*c1*c4=LG14p1*-2.0*c1*c4 +LG14p2*2.0*(c1-c4)*-2.0*c1*c4
					//d(LG24)*c2*c4=-LG24p*c2*c4
					//d(LG34)*c3*c4=-LG34p*c3*c4
					//=G0G4-G0G1
					//+LG12p*c1*c2
					//-LG13p*c1*c3
					//+LG14p1*-2.0*c1*c4 +LG14p2*2.0*(c1-c4)*-2.0*c1*c4
					//-LG24p*c2*c4
					//-LG34p*c3*c4
					//-LG12*c2 -LG13*c3 -LG14*c4
					//+LG14*c1 +LG24*c2 +LG34*c3
					//-LG123*c2*c3 -LG124*c2*c4+LG124*c1*c2 
					//+Gmg(omitted)+-RT*log(c1)-RT+RT*log(c4)+RT
					dgnig=beta*
						(G0G4-G0G1
						+LG12p*c1*c2
						-LG13p*c1*c3
						+LG14p1*-2.0*c1*c4 +LG14p2*2.0*(c1-c4)*-2.0*c1*c4
						-LG24p*c2*c4
						-LG34p*c3*c4
						-LG12*c2
						-LG13*c3
						-LG14*c4
						+LG14*c1
						+LG24*c2
						+LG34*c3
						-LG123*c2*c3
						-LG124*c2*c4
						+LG124*c1*c2
						+RT*(log(c4)-log(c1))
						);
					//----- ----- ----- ----- ----- -----
				
					/* Calculate the derivatives of free energy which
					   will be used later for the evaluation of 
					   derivative with respect to order parameter */
					//
					// Free energy G for Alpha and gamma, respectively
					//
					// G of alpha phase (modified original code)
					gcunia=beta*
						(G0A1*c1
						+G0A2*c2
						+G0A3*c3
						+G0A4*c4
						+LA12*c1*c2
						+LA13*c1*c3
						+LA14*c1*c4
						+LA23*c2*c3
						+LA24*c2*c4
						+LA34*c3*c4
						+LA123*c1*c2*c3
						+RT*(c1*log(c1)+c2*log(c2)+c3*log(c3)+c4*log(c4))
						);
					
					// G of gamma phase (modified original code)
					gcunig=beta*
						(G0G1*c1
						+G0G2*c2
						+G0G3*c3
						+G0G4*c4
						+LG12*c1*c2
						+LG13*c1*c3
						+LG14*c1*c4
						+LG23*c2*c3
						+LG24*c2*c4
						+LG34*c3*c4
						+LG123+c1*c2*c3
						+LG124*c1*c2*c4
						+RT*(c1*log(c1)+c2*log(c2)+c3*log(c3)+c4*log(c4))
						);
				}else{
					/* If Fe concentration is greater than 0.9995 or less than 0.0005 apply
					   a penalty term. This is a numerical trick to avoid the overflows
					   because of presence of log terms in free energy functional Eq.5.26 */
					
					if(c1>=0.9995){
						thrc=(c1-0.9995);
					}
					
					if(c1<=0.0005){
						thrc=(c1-0.0005);
					}
					
					//Calculate the derivatives of penalty term for each alloying elememnts
					gcuni=Coef*thrc*thrc;
					//----- -----
					gcunia=gcuni; //gcunia=Coef*thrc*thrc;
					gcunig=gcuni; //gcunig=Coef*thrc*thrc;
					//----- -----
					
					//Derivatives of free energy based on function h for each alloying elements
					dgcu=-2.0*Coef*thrc;
					//----- -----
					dgcua=dgcu; //dgcua=-2.0*Coef*thrc;
					dgcug=dgcu; //dgcug=-2.0*Coef*thrc;
					//
					dgmna=dgcu; //dgmna=-2.0*Coef*thrc;
					dgmng=dgcu; //dgmng=-2.0*Coef*thrc;
					//
					dgnia=dgcu; //dgnia=-2.0*Coef*thrc;
					dgnig=dgcu; //dgnig=-2.0*Coef*thrc;
					//----- -----
				}//end if
				
				//Derivative of free energy based on the function h for Cu,Mn,and,Ni concentration, respectively
				// d( [1-h(eta)]*[Galpha+Y*Vm*epsilon0(ci)^2] + h(eta)*Ggamma )/d(ci) and h(eta)=funch
				dgdcu[ii]=(1.0-funch)*(dgcua+delascu)+funch*dgcug;
				dgdmn[ii]=(1.0-funch)*(dgmna+delasmn)+funch*dgmng;
				dgdni[ii]=(1.0-funch)*(dgnia+delasni)+funch*dgnig;
				
				//Derivative of free energy with respsect to order parameter
				/* d( [1-h(eta)]*[Galpha+Y*Vm*epsilon0(ci)^2] + h(eta)*Ggamma + W*g(eta)^2)/d(eta)
				   h(eta)=eta*eta*(3-2*eta)=funch, g(eta)=eta*(1-eta)=funcg
				   d(h)/d(eta) = 2.0*eta*(3-2*eta) + eta*eta*-2.0 = 6.0*eta*(1.0-eta) = dhdor
				   d(g^2)/d(eta)=2.0*g*(1-2.0*eta) = 2.0*funcg*(1-2.0*eta) = funcg*(2.0-4.0*eta)
				   (gcunia+elaste) = [Galpha+Y*Vm*epsilon0(ci)^2]
				   gcunig = Ggamma, W = constw, eta=orp[ii]                                                  */
				//dgdor[ii][0]=(gcunia+elaste)*( 2.0*orp[ii]*(3.0-2.0*orp[ii]) + orp[ii]*orp[ii]*-2.0 )*-1.0
				//				  +gcunig*( 2.0*orp[ii]*(3.0-2.0*orp[ii]) + orp[ii]*orp[ii]*-2.0 )
				//				  +constw*( 2.0*orp[ii]*(1.0-orp[ii])*(1.0-2.0*orp[ii]) );
				//dhdor=6.0*orp[ii]*(1.0-orp[ii]);
				dgdor[ii]=(gcunia+elaste)*( -6.0*funcg ) // -(6.0*orp[ii]-6.0*orp[ii]*orp[ii])
								  +gcunig*(  6.0*funcg ) //  (6.0*orp[ii]-6.0*orp[ii]*orp[ii])
								  +constw*( funcg*(2.0-4.0*orp[ii]) );
			}//k
		}//j
	}//i
	return;
}