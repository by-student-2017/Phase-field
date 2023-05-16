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
	
	//Number of integration
	int ngaus2 = ngaus;
	if(nnode == 3){
		ngaus2 = 1;
	}
	mgaus = ngaus * ngaus2;
	
	for(int ielem=0;ielem<nelem;ielem++){
		//
		//Material Parameters and Elasticity Matrix
		mtype = matno[ielem];
		modps_2d(mtype,ntype,nstre,props,dmatx);
		poiss = props[mtype][2];
		
		//Nodal displacements
		for(int inode=0;inode<nnode;inode++){
			lnode = lnods[ielem][inode];
			for(int idofn=0;idofn<ndofn;idofn++){
				nposn = lnode*ndofn + idofn;
				eldis[idofn][inode] = asdis[nposn];
				elcod[idofn][inode] = coord[lnode][idofn];
			}
		}
		
		//Integrate Stresses
		kgasp = 0;
		for(int igaus=0;igaus<ngaus;igaus++){
			for(int jgaus=0;jgaus<ngaus;jgaus++){
				//
				kgasp = kgasp + 1;
				exisp = posgp[igaus];
				etasp = posgp[jgaus];
				
				if(ngaus2 == 1){
					etasp = posgp[ngaus+igaus];
				}
				
				sfr2_2d(exisp,stasp,nnode,shape,deriv);
				jacob2_2d(ielem,elcod,kgasp,shape,deriv,nnode,ndime,cartd,djacb,gpcod);
				bmats_2d(cartd,shape,inode,bmatx);
				
				//calculate the strains
				for(int istre=0;istre<nstre;istre++){
					strain[istre] = 0.0;
					for(int inode=0;inode<nnode;inode++){
						for(int idofn=0;idofn<ndofn;idofn++){
							ievab = lnode*ndofn + idofn;
							strain[istre] = strain[istre] + bmatx[istre][ievab] * eldis[idofn][inode];
						}
					}
				}
				
				//calculate stresses
				for(int istre=0;istre<nstre;istre++){
					stres[istre] = 0.0;
					for(int jstre=0;jstre<nstre;jstre++){
						stres[istre] = stres[istre] + dmatx[istre][jstre] * strain[jstre];
					}
				}
				
				if(ntype == 1){
					stres[3] = 0.0;
				}
				if(ntype ==2){
					stres[3] = poiss*(stres[0]+stres[1]);
				}
				
				for(int istre=0;istre<nstre+1;istre++){
					elem_stres[ielem][kgasp][istre] = stres[istre];
				}
				//
			}//end for(igaus
		}//end for(jgaus
		//
	}//end for(ielem
	
	return;
}