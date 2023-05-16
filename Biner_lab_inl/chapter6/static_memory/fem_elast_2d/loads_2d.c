/* This function integrates the prescribed distributed external loads at
   the element edges and also the prescribed point loads at
   the nodes into the global force vector. */

#include <stdio.h> //printf()

/* Variable and array list
  npoin: Total number of nodes in the solution
  nelem: Total number of element in the solution
  ndofn: Number of DOFs per node
  ngaus: The order of numerical integration
  ndime: Number of global Cartesian coordinates
  eload[nelem][nevab]: Element force vector (nevab = nnode * ndofn).
  lnods[nelem][nnode]: Element nodal connectivity list
  coord[npoin][ndime]: Cartesian coordinates of nodes
  posgp[ngaus]: Position of sampling points for numerical integration
  weigp[ngaus]: Weight of sampling points for numerical integration
  gforce[ntotv]: Global force vecotr (ntotv = npoin * ndofn).
*/

void sfr2_2d();

void loads_2d(int npoin, int nelem, int ndofn, 
	int nnode, int ngaus, int ndime,
	double *posgp, double *weigp, double *lnods, double *coord,
	double *gforce, FILE *in, FILE *out){
	
	//Initilize global force vector (rhs) & element loads
	nevab = nnode * ndofn;
	ntotv = npoin * ndofn;
	
	for(int itotv=0;itotv<ntotv;itotv++){
		gforce[itotv] = 0.0;
	}
	
	for(int ielem=0;ielem<nelem;ielem++){
		for(int ievab=0;ievab<nevab;ievab++){
			eload[ielem][ievab] = 0.0;
		}
	}
	
	//Loading types
	fscanf(in,"%5d %5d",&iplod,&nedge);
	
	//Point forces
	if(iplod != 0){
		fscanf(in,"%5d",&nplod);
		for(int jplod=0;jplod<nplod;jplod++){
			//
			fscanf(in,"%5d",&lodpt);
			for(int idofn=0;idofn<ndofn;idofn++){
				fscanf(in,"%lf",&point[idofn]);
			}
			//
			for(int ielem=0;ielem<nelem;ielem++){
				for(int inode=0;inode<nnode;inode++){
					nloca = lnods[ielem][inode];
					if(lodpt == nloca){
						for(int idofn=0;idofn<ndofn;idofn++){
							nposi = (inode * ndofn) + idofn;
							eload[ielem][nposi] = point[idofn];
						}
					}
				}
			}
			//
		}
	}//end if(iplod
	
	//Distributed Forces
	if(nedge != 0){
		fprintf(out, "Number of loaded edges: %5d \n",nedge);
		fprintf(out, "List of loaded edges and applied loads:\n");
		
		for(int iedge=0;iedge<nedge;iedge++){
			//
			nodeg = 3;
			if(nnode != 8){
				nodeg = 2;
			}
			
			//Read loads
			fscanf(in,"%5d",&neass);
				for(int iodeg=0;iodeg<nodeg;iodeg++){
				fscanf(in,"%5d",&noprs[iodeg]);
			}
			
			for(int iodeg=0;iodeg<nodeg;iodeg++){
				for(int idofn=0;idofn<ndofn;idofn++){
					fscanf(in,"%lf",&press[iodeg][idofn]);
				}
			}
			
			//Print
			fprintf(out, "\n");
			fprintf(out, "%5d \n",neass);
			for(int iodeg=0;iodeg<nodeg;iodeg++){
				fprintf(out, "%5d \n",noprs[iodeg]);
			}
			fprintf(out, "\n");
			for(int iodeg=0;iodeg<nodeg;iodeg++){
				for(int idofn=0;idofn<ndofn;idofn++){
					fprintf(out, "%14.6e \n",press[iodeg][idofn]);
				}
			}
			
			//End of reading
			
			//Integrated along the edges
			
			etasp = -1.0;
			for(int iodeg=0;iodeg<nodeg;iodeg++){
				lnode = noprs[iodeg];
				for(int idime=0;idime<ndime;idime++){
					elcod[idime][iodeg] = coord[lnode][idime];
				}
			}
			
			for(int igaus=0;igasu<ngaus;igaus++){
				exisp = posgp[igaus];
				sfr2_2d(exisp,etasp,nnode,shape,deriv);
				
				for(int idofn=0;idofn<ndofn;idofn++){
					pgash[idofn] = 0.0;
					dgash[idofn] = 0.0;
					for(int iodeg=0;iodeg<nodeg;iodeg++){
						pgash[idofn] = pgash[idofn]
									 + press[iodeg][idofn] * shape[iodeg];
						dgash[idofn] = dgash[idofn]
									 + elcod[idofn][iodeg] * deriv[0][idoeg];
					}
				}
				
				dvolu = weigp[igasu];
				pxcom = dgash[0] * pgash[1] - dgash[1] * pgash[0];
				pycom = dgash[0] * pgash[0] + dgash[1] * pgash[1];
				
				for(int inode=0;inode<node;inode++){
					nloca = lnods[neass][inode];
					if(nloca == noprs[0]){
						kount = 0;
						for(int knode=inode;knode<jnode;knode++){
							kount = kount + 1;
							ngash = (knode * ndofn) + 0;
							mgash = (knode * ndofn) + 1;
							if(knode > nnode){
								ngash = 0;
								mgash = 1;
							}
							eload[neass][ngash] = eload[neass][nagsh]
												+ pxcom * dvolu * shape[kount];
							eload[neass][mgash] = eload[neass][magsh]
												+ pycom * dvolu * shape[kount];
						}
					}
				}
			}//end for(igauss
		}//end for(iedge
		//
	}//end for(nedge
	
	//Print Nodal Forces
	fprintf(out, "\n");
	fprintf(out, "Nodal forces for elements: \n");
	for(int ielem=0;ielem<nelem;ielem++){
		fprintf(out, "\n");
		fprintf(out, "Element No: %5d \n",ielem);
		for(int ievab=0;ievab<nevab;ievab++){
			fprintf(out, "%14.6e ",eload[ielem][ievab]);
			if( (nnode == 8) && (ievab == nevab/2) ){
				fprintf(out, "\n");
			}
		}
		fprintf(out, "\n");
	}
	
	//Generate global force vector
	for(int ielem=0;ielem<nelem;ielem++){
		for(int inode=0;inode<nnode;inode++){
			lnode = lnods[ielem][inode];
			for(int idofn=0;idofn<ndofn;idofn++){
				itotv = (lnode * ndofn) + idofn;
				ievab = (inode * ndofn) + idofn;
				gforce[itotv] = gforce[itotv] + eload[ielem][ievab];
			}
		}
	}
	
	return;
}