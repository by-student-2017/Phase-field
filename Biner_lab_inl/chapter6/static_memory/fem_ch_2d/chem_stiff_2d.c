/* This function forms the global stiffness and global rhs, load, vector for solution of
   Cahn-Hilliard equation for the modified Newton-Raphson solultion algorithm. */

/* Variable and array list
  npoin: Total number of nodes in the solution (int)
  nelem: Total number of elements in the solution (int)
  nnode: Number of nodes per element (int)
  nstre: Number of stress compoents (int)
  ndime: Number of Cartesian coordinate dimension (int)
  ndofn: Number of DOF per node (int)
  ngaus: The order of numerical integration (int)
  ntype: Solution type (ntype=1 for plane-stress and,
    ntype=2 for plane-strain) (int)
  mobil: Mobility coefficient (double)
  grcoef: Gradient energy coefficient (double)
  dtime: Time increment used in time integration (double)
  istep: Time increment step number (int)
  iter: Newton-Raphson iteration number (int)
  con[ntotv]: Nodal values, which includes c and mu values,
    int ntotv = npoin * ndofn. (double)
  con_old[ntotv]: Nodal values from previous time step (double)
  lnods[nelem][nnode]: Element nodal connectivity list (int)
  coord[npoin][ndime]: Cartesian coordinates of nodes (double)
  posgp[ngaus]: Position of sampling points (double)
  weigp[ngaus]: Weighting factors at the sampling points (double)
  gstif[ntotv][ntotv]: Global stiffness matrix (double)
  gforce[ntotv]: Global rhs, load, vector (double)
*/

void sfr2_2d();
void jacob2_2d();
void free_energy_fem_2d();

void chem_stiff_2d(int npoin, int nelem, int nnode, int nstre, int ndime,
	int ndofn, int ngaus, int ntype,
	int *lnods, double *coord,
	double mobil, double grcoef,
	double *con, double *con_old,
	double dtime,
	double *posgp, double *weigp,
	int istep, int iter,
	double *gstif, double *gforce){
	
	//----- ----- ----- ----- ----- -----
	// Global and local variables
	//----- ----- ----- ----- ----- -----
	int ntotv = npoin * ndofn;
	int nevab = nnode * ndofn;
	//----- 
	int ngaus2 = ngaus;
	
	if(nnode == 3){
		ngaus2 = 1;
	}
	
	for(int i=0;i<ntotv;i++){
		gforce[i] = 0.0;
	}
	
	//----- ----- ----- ----- ----- -----
	// initialize elements stiffness & rhs
	//----- ----- ----- ----- ----- -----
	// stiffness matrices
	double kcc[nnode][nnode];
	double kcm[nnode][nnode];
	double kmm[nnode][nnode];
	double kmc1[nnode][nnode];
	double kmc[nnode][nnode];
	if(iter == 1){
		for(int inode=0;inode<nnode;inode++){
			for(int jnode=0;jnode<nnode;jnode++){
				kcc[inode][jnode] = 0.0;
				kcm[inode][jnode] = 0.0;
				kmm[inode][jnode] = 0.0;
				kmc1[inode][jnode] = 0.0;
				kmc[inode][jnode] = 0.0;
			}
		}
	}//end if
	
	//rhs
	double eload1[nnode];
	double eload2[nnode];
	for(int inode=0;inode<nnode;inode++){
		eload1[inode] = 0.0;
		eload2[inode] = 0.0;
	}
	double eload[nevab];
	for(int ievab=0;ievab<nevab;ievab++){
		eload[ievab] = 0.0;
	}
	
	//----- ----- ----- ----- ----- -----
	// elemental values
	//----- ----- ----- ----- ----- -----
	double cv[nnode];
	double cm[nnode];
	double cv_old[nnode];
	for(int inode=0;inode<nnode;inode++){
		lnode = lnods[ielem][inode];
		cv[inode] = con[lnode];
		cm[inode] = con[npoin + lnode];
		cv_old[inode] = con_old[lnode];
	}
	
	//coords of the element nodes
	for(int inode=0;inode<nnode;inode++){
		lnode = lnods[ielem][inode];
		for(int idime=0;idime<ndime;idime++){
			elcod[idime][inode] = coord[lnode][idime];
		}
	}
	
	//----- ----- ----- ----- ----- -----
	// integrate elements stiffness and rhs
	//----- ----- ----- ----- ----- -----
	int kgasp = 0;
	
	for(int igaus=0;igasu<ngasu;igaus++){
		exisp = posgp[igaus];
		for(int jgaus=0;jgasu<ngasu;jgaus++){
			etasp = posgp[jgaus];
			if(nnode == 3){
				etasp = posgp[ngaus + igaus];
			}
			
			kgasp = kgasp + 1;
			sfr2_2d(exisp,etasp,nnode,shape,deriv);
			jacob2_2d(ielem,elcod,kgasp,shape,
				deriv,nnode,ndime,cartd,djacb,gpcod);
			
			dvolu = djacb * weigp[igaus] * weigp[jgaus];
			
			
		}
	}
	
	return;
}