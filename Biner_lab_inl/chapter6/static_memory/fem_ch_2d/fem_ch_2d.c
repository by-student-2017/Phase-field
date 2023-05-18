/* Finite Element Method (FEM) code for solving Cahn-Hilliard equation */

/* This is the main program to solve Chan-Hilliard phase-field equation with FEM algorithm.  */

#include <stdio.h> //printf()
#include <stdlib.h> //rand()
#include <math.h> //mod() and -lm
#include <time.h>

void input_fem_pf_2d();
void periodic_boundary_2d();
void gauss_2d();
void cart_deriv_2d();
void init_micro_ch_fem_2d();
void chem_stiff_2d();
void apply_periodic_bc_2d();
void recover_slave_dof_2d();
void write_vtk_fem_2d();

int main(){
	// Get initial wall clock time beginning of the execution
	clock_t start, end;
	double compute_time;
	//get initial wall time
	start = clock();
	
	//open input files
	FILE *in=fopen("mesh3n_40.inp","r");
	
	//open output files
	FILE *out=fopen("result_1.out","w");
	
	//int isolve=1;
	
	//time integration parameters
	int nstep = 5000; //Number of time increment steps
	int nprnt = 25;   //Print frequency to output the results to file
	double dtime = 2.0e-2; //Time increment for numerical integration
	double toler = 5.0e-5; //Tolerance value for iterative solution
	int miter = 10; //Number of maximum iterations
	
	//Material specific parameters
	double conc0 = 0.40; //Alloy concentration
	double mobil = 1.0;  //Mobility parameter
	double grcoef = 0.5; //Gradient energy coefficient
	
	//----- ----- ----- ----- ----- -----
	// input data
	//----- ----- ----- ----- ----- -----
	// Input FEM mesh and control parameters
	input_fem_pf_2d(npoin,nelem,ntype,nnode,ndofn,ndime,ngaus,nstre,nmats,nprop,lnods,matno,coord);
	
	// Determine the master and slave nodes for periodic boundary condition
	periodic_boundary_2d(npoin,coord,ncouontm,ncounts,master,slave);
	
	//Parameters for numercial integration of elements
	gauss_2d(ngaus,nnode,posgp,weigp);
	
	int ntotv = npoin * ndofn; //Total number of variables in the solution
	int nevab = nnode * ndofn; //Total number of variables per element
	
	//----- ----- ----- ----- ----- -----
	// Prepare microstructure
	//----- ----- ----- ----- ----- -----
	// Modulate the microstructure and initialize the vector containing nodal variables
	init_micro_ch_fem_2d(npoin,ndofn,conc0,con);
	
	//----- ----- ----- ----- ----- -----
	// Evolve
	//----- ----- ----- ----- ----- -----
	// Evolve microstructure
	for(int istep=0;istep<nstep;istep++){
		
		// Assign nodal variables ffrom previous time step
		for(int i=0;i<ntotv;i++){
			con_old[i] = con[i];
		}
		
		//----- ----- ----- ----- ----- -----
		// Newton interation
		//----- ----- ----- ----- ----- -----
		// Modified Newton-Raphson iterations
		for(int iter=0;iter<miter;iter++){
			
			// If iter=1, initalize the global stiffness matrix
			if(iter==0){
				sparse_2d(ntotv,ntotv,gstif);
			}
			
			/* Form the global stiffness matrix for only in
			   first iteration and form rhs, load, vector in
			   every following iterations. */
			chem_stiff_2d(npoin,nelem,nnode,nstre,ndime,gstif,gforce,
				ndofn,ngaus,ntype,lnods,coord,
				mobil,grcoef,con,con_old,dtime,
				posgp,weigp,istep,iter,gstif);
			
			//----- -----
			// Rearrange gstif & gforce for PBC
			//----- -----
			/* Modify the global stiffness matrix and rhs,
			   load, vector for periodic boundaries. */
			apply_periodic_bc(ncountm,ncounts,master,slave,
				ndofn,npoin,gstif,gforce,iter);
			
			//----- -----
			// solve squation and update
			//----- -----
			// Solve the resulting linear equations
			//----- ----- ----- ----- ----- -----
			// Solve displacements (asdis = gstif\gforce)
			//  left division: inverse(gstif)*gforce
			// http://www.eccse.kobe-u.ac.jp/assets/files/2020/Kobe_HPC_Spring_terao.pdf
			// https://ist.ksc.kwansei.ac.jp/~nishitani/Lectures/2005/NumRecipeEx/Lapack.pdf
			//----- ----- ----- ----- ----- -----
			// Solve the resulting system of equation for the nodal displacements
			//
			// Transpose
			double gstif_t[ntotv][ntotv];
			for(int x=0;x<ntotv;x++){
				for(int y=0;y<ntotv;y++){
					gstif_t[y][x] = gstif[x][y];
				}
			}
			//
			double asdis[ntotv];
			int n=ntotv, nrhs=1, lda=ntotv, ldb=ntotv, info=0;
			dgesv_(&n, &nrhs, gstif_t, &lda, asdis, gforce, &ldb, &info); //x=A\B for Ax=B, C language -> Fortran library
			
			//----- ----- ----- ----- ----- -----
			// Recover slave node values
			//----- ----- ----- ----- ----- -----
			// Recover the nodal values for the slave nodes in the periodic boundaries.
			recover_slave_dof_2d(asdis,ncountm,ncounts,master,slave,npoin,ndofn);
			
			//update concentration field
			// Update the nodal ariable vector
			for(int i=0;i<ntotv;i++){
				con[i] = con[i] + asdis[i];
			}
			
			//for small deviations
			// For small derivations from max and min values, reset the limits.
			for(int ipoin=0;ipoin<npoin;ipoin++){
				if(con[ipoin] >= 0.9999){
					con[ipoin] = 0.9999;
				}
				if(con[ipoin] <= 0.0001){
					con[ipoin] = 0.0001;
				}
			}
			
			//check norm for convergence
			/* Check convergence. If converged solution is reached,
			   exit from Newton-Raphson iteration. */
			//normF = norm(gforce,2);
			double normF = 0.0;
			for(int i=0;i<ntotv;i++){
				normF = normF + gforce[i]*gforce[i];
			}
			
			if(normF <= toler){
				break;
			}
			//
		}//end for(iter
		
		//print out
		// If print frequency is reached, output the results to file
		if(fmod(istep,nprint)==0){
			printf("done step: %5d \n",istep);
			
			// Write the results in vtk format for contour plots to be veiwed by Paraview.
			write_vtk_fem_2d(npoin,nelem,nnode,lnods,coord,istep,con);
		}
		//
	}//end for(istep
	
	//Display total execution time for the solution
	// Compute the total execution time and print
	end = clock();
	compute_time = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("Compute Time: %lf \n", compute_time);
	
	close(in);
	close(out);
}
