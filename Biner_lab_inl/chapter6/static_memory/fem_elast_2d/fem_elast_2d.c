/* Finite Element Method (FEM) code for linear elasticity */

/* This FEM program solves the linear elasticity problems in
   two dimensions by using three-, four-, and eigth-node
   isoparametric elements. */

#include <stdio.h> //printf()
#include <stdlib.h> //rand()
#include <math.h> //mod() and -lm
#include <time.h>

void input_fem_elast_2d();
void gauss_2d();
void stiffness_2d();
void loads_2d();
void boundary_cond_2d();
void stress_2d();
void output_fem_2d();

int main(){
	// Get initial wall clock time beginning of the execution
	clock_t start, end;
	double compute_time;
	//get initial wall time
	start = clock();
	
	//open input files
	FILE *in=fopen("mesh_1.inp","r");
	
	//open output files
	FILE *out=fopen("result_1.out","w");
	
	//----- ----- ----- ----- ----- -----
	// Difinition (Variable and array list)
	//----- ----- ----- ----- ----- -----
	int npoin; //Total number of nodes in the solution
	int nelem; //Total number of elements
	int nvfix; //Total number of constrained nodes
	int ntype; //Solution type (ntyep=1, plane-stress and ntype=2 plane-strain)
	int nnode; //Number of nodes per element
	int ndofn; //Number of degree of freedom (DOF) per node
	int ndime; //Number of spatial dimensions
	int ngaus; //The order of numerical integration
	int nmats; //Total number of different materials in the solution
	int nstre; //Number of stress components
	int nprop; //Number of material properties
	//
	int matno[nelem]; //Material types for the elements
	int nofix[nvfix]; //Node numbers at which one or more DOFs are constrained
	//
	int    lnods[nelem][nnode]; //Element nodal connectivity list
	double coord[npoin][ndime]; //Cartesian coordinates of each node
	int    iffix[nvfix][ndofn]; //List of constrained DOFs
	double fixed[nvfix][ndofn]; //Prescribed value of any constrained DOFs
	double props[nmats][nprop]; //For each different material, the properties of that material
	
	//----- ----- ----- ----- ----- -----
	// Input FEM data
	//----- ----- ----- ----- ----- -----
	input_fem_elast_2d(npoin,nelem,nvfix,ntype,nnode,ndofn,ndime,ngaus,nmats,nstre,nprop,
		matno,nofix,lnods,coord,iffix,fixed,props,
		in,out); // Input solution variables and FEM mesh
	int ntotv = npoin * ndofn;
	double posgp[ngaus];
	double weigp[ngaus];
	gauss_2d(ngaus,nnode,posgp,weigp); // Get the values for the numerical integration
	
	//----- ----- ----- ----- ----- -----
	// From global stifnes matrix (lhs)
	//----- ----- ----- ----- ----- -----
	// From the global stifness matrix
	double gstif[ntotv][ntotv];
	stiffness_2d(npoin,nelem,nnode,nstre,ndime,ndofn,ngaus,
		ntype,lnods,matno,coord,props,posgp,weigp,gstif);
	
	//----- ----- ----- ----- ----- -----
	// Force vector & Boundary Conditions
	//----- ----- ----- ----- ----- -----
	double gforce[ntotv]; // ntotv = npoin * ndofn;
	loads_2d(npoin,nelem,ndofn,nnode,ngaus,
		ndime,posgp,weigp,lnods,coord,gforce); // Form the global force vector from the loading information
	boundary_cond_2d(npoin,nvfix,nofix,iffix,fixed,ndofn,gstif,gforce); // Apply boundary conditions
	
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
	// Calculate stresses
	//----- ----- ----- ----- ----- -----
	// Calculate stress values
	stress_2d(asdis,nelem,npoin,nnode,ngaus,
		nstre,props,ntype,ndofn,ndime,lnods,matno,coord,posgp,weigp);
	
	//----- ----- ----- ----- ----- -----
	// Output results
	//----- ----- ----- ----- ----- -----
	// Output displacements and stress values to files
	double elem_stres[nelem][ngaus][nstre];
	output_fem_2d(npoin,nelem,nnode,lnods,coord,ndofn,ngaus,nstre,asdis,elem_stres);
	
	//end of execution
	printf("Done \n");
	
	//Display total execution time for the solution
	end = clock();
	compute_time = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("Compute Time: %lf \n", compute_time);
	
	close(in);
	close(out);
}
