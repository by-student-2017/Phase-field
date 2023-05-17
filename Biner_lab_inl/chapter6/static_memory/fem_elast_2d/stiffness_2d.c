/* This function forms the stiffness matrices at
   the element level first, then, assembles them into
   the global stiffness matrix. */

/* Variable and array list
  npoin: Total number of nodes in the solution (int)
  nelem: Total number of elements in the solution (int)
  nnode: Number of nodes per element (int)
  nstre: Number of stress components (int)
  ndime: Number of global Cartesian components (int)
  ndofn: Number of degrees of feedom per node (int)
  ngaus: The order of Gaussian integration (int)
  ntype: solution type, ntype=1 for plane-stress,
    ntype=2 for plane-strain (int)
  matno[nelem]: Material types for elements (int)
  posgp[ngaus]: Position of sampling points for numerical integration. (double)
  weigp[ngaus]: Weights of sampling points for numerical integration. (double)
  lnods[nelem][nnode]: Element nodal connectivity list (int)
  coord[npoin][ndime]: Cartesian coordinates of nodes (double)
  props[matno][ndime]: For each different material,
    the properties of that material (double)
  gstif[ntotv][ntotv]: Global stiffness matrix (ntotv = npoin * ndofn). (double)
*/

void sfr2_2d();
void jacob2_2d();
void modps_2d();
void bmats_2d();
void dbe_2d();

void stiffness_2d(int npoin, int nelem, int nnode,
	int nstre, int ndime, int ndofn, int ngaus, int ntype,
	int *lnods, int *matno, double *coord,
	double *props, double *posgp, double *weigp, double *gstif){
	
	//----- ----- ----- ----- ----- ----- 
	//Initialize the global stiffness
	
	// Calculate total number of element variables
	int nevab = nnode * ndofn;
	
	/* Change the order of numerical integration parameter for
	   three-node isoparametric element. */
	int ngaus2 = ngaus;
	if(nnode == 3){
		ngaus2 = 1;
	}
	
	/* Calculate the total number of variables in
	   the solution and initialize the global stiffness matrix. */
	int ntotv = npoin * ndofn;
	for(int i=0;i<ntotv;i++){
		for(int j=0;j<ntotv;j++){
			gstif[i*ntotv+j] = 0.0;
		}
	}
	//----- ----- ----- ----- ----- ----- 
	
	double estif[nevab][nevab];
	int mtype;
	double dmatx[nstre][nstre];
	int lnode;
	double elcod[ndime][nnode];
	int kgasp;
	double exisp, etasp;
	double shape[nnode];
	double deriv[ndime][nnode];
	double cartd[ndime][nnode];
	double djacb;
	double gpcod[ndime][ngaus*ngaus];
	double bmatx[nstre][nevab];
	double dbmat[nstre][nevab];
	double dvolu;
	int itotv, jtotv, ievab, jevab, knode;
	//----- ----- ----- ----- ----- ----- 
	//Element stiffness & loads
	for(int ielem=0;ielem<nelem;ielem++){
		
		//Initialize element stiffness
		// Initialize element stiffness matrix
		for(int ievab=0;ievab<nevab;ievab++){
			for(int jevab=0;jevab<nevab;jevab++){
				estif[ievab][jevab] = 0.0;
			}
		}
		
		//Form elasticity matrix
		/* Depending upon the material type of
		   current element form the elasticity matrix Eqs.6.39 or 6.40. */
		mtype = matno[ielem];
		modps_2d(mtype,ntype,nstre,props,dmatx);
		
		//Coordinates of element nodes
		/* define coordinates of the nodes of 
		   the current element (ielem) in a local array. */
		for(int inode=0;inode<nnode;inode++){
			lnode = lnods[ielem*nnode+inode];
			for(int idime=0;idime<ndime;idime++){
				elcod[idime][inode] = coord[lnode*ndime+idime];
			}
		}
		
		//Integrate the element stiffness
		// Numerical integration loop
		kgasp = 0; //Counter for integration points
		// Determine the positions of integration points
		for(int igaus=0;igaus<ngaus;igaus++){
			exisp = posgp[igaus];
			for(int jgaus=0;jgaus<ngaus;jgaus++){
				etasp = posgp[jgaus];
				if(nnode == 3){
					etasp = posgp[ngaus + igaus];
				}
				
				// Increment integration point counter
				kgasp = kgasp + 1;
				
				/* Calculate the thape function and 
				   their derivatives for the current integration point. */
				sfr2_2d(exisp,etasp,nnode,shape,deriv);
				
				/* Calculate the Cartesian derivatives of the shape functions,
				   the determinant of the Jacobian, and the cartesian coordinates of
				   the integration points (Eq.6.6). */
				jacob2_2d(ielem,elcod,kgasp,shape,deriv,nnode,ndime,cartd,djacb,gpcod);
				
				// Calculate the strain matrix (Eq.6.37)
				bmats_2d(cartd,shape,nnode,bmatx);  //inode ?
				
				// Multiply the elasticity matrix with strain matrix
				dbe_2d(nevab,nstre,bmatx,dmatx,dbmat);
				
				// Calculate the area of the element (Eq.6.7)
				dvolu = djacb * weigp[igaus] * weigp[jgaus];
				
				/* Correct the area of the element, if ielem is
				   a three-node isoparametric element. */
				if(nnode == 3){
					dvolu = djacb * weigp[igaus];
				}
				
				//From element stiffness
				// From the element stiffness matrix (Eq.6.7)
				for(int ievab=0;ievab<nevab;ievab++){
					for(int jevab=0;jevab<nevab;jevab++){
						for(int istre=0;istre<nstre;istre++){
							estif[ievab][jevab] = estif[ievab][jevab]
												+ bmatx[istre][ievab] * dbmat[istre][jevab] * dvolu;
						}
					}
				}
			}//end for(jgaus
		}//end for(igasu
		
		// From Global stiffness matrix
		/* Place the components of the element stiffness matrix into
		   global stiffness matrix based on the information provided by
		   the elemental connectivity list, lnods. */
		for(int inode=0;inode<nnode;inode++){
			lnode = lnods[ielem*nnode+inode];
			for(int idofn=0;idofn<ndofn;idofn++){
				itotv = (lnode * ndofn) + idofn;
				ievab = (inode * ndofn) + idofn;
				//
				for(int jnode=0;jnode<nnode;jnode++){
					knode = lnods[ielem*nnode+jnode];
					for(int jdofn=0;jdofn<ndofn;jdofn++){
						jtotv = (knode * ndofn) + jdofn;
						jevab = (jnode * ndofn) + jdofn;
						//
						gstif[itotv*ntotv+jtotv] = gstif[itotv*ntotv+jtotv]
												 + estif[ievab][jevab];
					}//end for(jdofn
				}//end for(jnode
				//
			}//end for(idofn
		}//end for(inode
		
	}//end for(ielem
	//----- ----- ----- ----- ----- ----- 
	
	return;
}