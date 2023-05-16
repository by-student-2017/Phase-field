/* This function rearranges the stiffness matrix and
   force vector for prescribed boundary conditions */

#include <stdio.h> //printf()

/* Variable and array list
  npoin: Total number of nodes in the solution
  nvfix: Total number of nodes with prescribed boundary conditions
  ndofn: Number of DOF per node
  nofix[nvfix]: Node numbers of which onw or more DOF are constrained
  fixed[nvfix]: The values of prescribed DOFs
  gforce[ntotv]: Global force vector, ntotv = npoin * ndofn
  iffix[nvfix][ndofn]: List of constrained DOFs
  gstif[ntotv][ntotv]: Global stiffness matrix
*/

void boundary_cond_2d(int npoin, int nvfix,
	int *nofix, double *fixed,
	int *ndofn, double *gstif, double *gforce,
	double *gstif, double *gforce){
	
	ntotv = npoin * ndofn;
	
	for(int ivfix=0;ivfix<nvfix;ivfix++){
		lnode = nofix[ivfix];
		for(int idofn=0;idofn<ndofn;idofn++){
			if(iffix[ivfix][idofn] == 1){
				itotv = lnode*ndofn + idofn;
				
				for(int jtotv=0;jtotv<ntotv;jtotv++){
					gstif[itotv][jtotv] = 0.0;
				}
				
				gstif[itotv][itotv] = 1.0;
				gforce[itotv] = fixed[ivfix][idofn];
			}//end if
		}//end for(idofn
	}//end for(ivfix
	
	return;
}