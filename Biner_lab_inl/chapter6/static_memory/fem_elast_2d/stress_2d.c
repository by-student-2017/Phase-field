/* This function calculates the strain and stress values at
   the integration points of the elements resulting from
   the displacements */

#include <stdio.h> //printf()

/* Variable and array list
  nelem: Total number of elements in the solution
  npoin: total number of the nodes in the solution
  nnode: Number of nodes per element
  ngaus: The order of numerical integration
  nstre: Number of stress components
  ntype: Solution type (ntype=1 for plane-stress,
    ntype=2 for plane-strain)
  ndofn: Number of DOF per node
  ndime: Number of global Cartesian coordinate components
  asdis[ntotv]: Global displacement vector (ntotv = npoin * ndofn)
  matno[nelem]: Material types for the elements
  posgp[ngaus]: Position of sampling points for numerical integration
  weigp[ngaus]: Weights of sampling points for numerical integration
  lnode[nelem][nnode]: Element nodal connectivity list
  props[nmats][nprop]: For each different material the properties of that material
  coord[npoin][dnime]: Cartesian coordinates of nodes
  elem_stres[nelem][mgaus][nstre]: stress componennts at element gauss points
*/

void modps_2d();
void sfr2_2d();
void jacob2_2d();
void bmats_2d();

void stress_2d(int asdis, int nelem, int npoin, int nnode,
	int ngaus, int nstre,
	double *props, int ntype, int ndofn, int ndime,
	int lnode, int *matno, 
	double *coord, double *posgp, double *weigp,
	double elem_stres){
	
	
	
	return;
}