/* This function prints out the solution results to file in
   tabulated form and in vtk file format for viewing the results by
   using Paraview */

#include <stdio.h> //printf()

/* Variable and array list
  npoin: Total number of nodes in the solution (int)
  nelem: total number of elements in the solution (int)
  nnode: Number of nodes per element (int)
  ndofn: Number of dgree of freedom per node (int)
  ngaus: The order of numerical integration (int)
  nstre: Number of stress components (int)
  asdis[ntotv]: Global displacement vector, ntotv = npoin * ndofn. (double)
  lnods[nelem][nnode]: Element connectivity list (int)
  coord[npoin][ndime]: Cartesian coordinates of the nodes (double)
  elem_stres[nelem][mgaus][nstre]: Element stress values at integration points (double)
*/

void write_vtk_fem_2d();

void output_fem_2d(int npoin, int nelem, int nnode,
	int *lnods, double *coord, int ndofn, 
	int ngaus, int nstre,
	double *asdis, double elem_stres,
	FILE *out){
	
	fprintf(out, "\n");
	fprintf(out, "******************\n");
	fprintf(out, "* outputs \n");
	fprintf(out, "******************\n");
	
	//Print Nodal Displacements
	// Write the displacement values tabulated format to the output file
	fprintf(out, "\n");
	fprintf(out, "Displacements:\n");
	fprintf(out, "Node Number X-Disp Y-Disp \n");
	
	for(int ipoin=0;ipoin<npoin;ipoin++){
		fprintf(out, "%5d \n",ipoin);
		for(int idofn=0;idofn<ndofn;idofn++){
			itotv = ipoin*ndofn + idofn;
			fprintf(out, "%14.6e \n",asdis[itotv]);
		}
		fprintf(out, "\n");
	}
	
	//Print Stress values
	// Write the element stress values in tabulated format to the output file
	fprintf(out, "\n");
	fprintf(out, "Stress at elements \n");
	fprintf(out, "gauss-point sigma-11 sigma-22 sigma-12 sigma-33 \n");
	
	//Number of integration
	ngaus2 = ngaus;
	if(nnode == 3){
		ngaus2 = 1;
	}
	
	for(int ielem=0;ielem<nelem;ielem++){
		fprintf(out, "\n");
		fprintf(out, "Element No: %5d \n",ielem);
		
		kgasp = 0;
		for(int igasu=0;igasu<ngasu;igasu++){
			for(int jgaus=0;jgasu<ngasu;jgasu++){
				kgasp = kgasp + 1;
				fprintf(out, "%d %14.6e %14.6e %14.6e %14.6e \n",
					kgasp,elem_stres[ielem][kgasp][0],
					elem_stres[ielem][kgasp][1],elem_stres[ielem][kgasp][2],elem_stres[ielem][kgasp][3],
					);
			}
		}
	}//end for(ielem
	
	//prepare results for graphical output (vtk file format)
	// deformed mesh
	/* Add displacement component values, amplified with the value of the facto,
	   to the nodal coordinates, if viewing of deformed mesh is desired. */
	facto = 3.0;
	for(int ipoin=0;ipoin<npoin;ipoin++){
		for(int idofn=0;idofn<ndofn;idofn++){
			itotv = ipoin*ndofn + idofn;
			disp_cord[ipoin][idofn] = coord[ipoin][idofn] + facto * asdis[itotv];
		}
	}
	
	//Extraplot stresses from integration points and average over entire mesh
	// number of connection of nodes
	/* Find out for each node, their number of connectivity. */
	for(int ipoin=0;ipoin<npoin;ipoin++){
		node_con[ipoin] = 0;
		for(int ielem=0;ielem<nelem;ielem++){
			for(int inode=0;inode<nnode;inode++){
				lnode = lnods[ielem][inode];
				if(lnode == ipoin){
					node_con[ipoin] = node_con[ipoin] + 1;
				}
			}
		}
	}//end for(ipoin
	
	//initialize nodal_stress
	/* Accumulate the stress values at the nodes, 
	   from element integration points based on the element connectivity list. */
	for(int ipoin=0;ipoin<npoin;ipoin++){
		for(int istre=0;istre<nstre;istre++){
			node_stres[ipoin][istre] = 0.0;
		}
	}
	
	for(int ielem=0;ielem<nelem;ielem++){
		//
		for(int istre=0;istre<nstre;istre++){
			ave_stres[istre] = 0.0;
		}
		
		kgasp = 0.0;
		for(int igaus=0;igasu<ngasu;igasu++){
			for(int jgaus=0;jgasu<ngasu;jgasu++){
				kgasp = kgasp + 1;
				for(int istre=0;istre<nstre;istre++){
					ave_stres[istre] = ave_stres[istre] + elem_stres[ielem][kgasp][istre];
				}
			}
		}
		
		for(int inode=0;inode<nnode;inode++){
			lnode = lnods[ielem][inode];
			for(int istre=0;istre<nstre;istre++){
				node_stre[lnode][istre] = node_stres[lnode][istre] + ave_stres[istre] / kgasp;
			}
		}
		//
	}//end (ielem
	
	// Average the nodal stress values based on their number of connections.
	for(int ipoin=0;ipoin<npoin;ipoin++){
		for(int istre=0;istre<nstre;istre++){
			node_stres[ipoin][istre] = node_stres[ipoin][istre] / node_con[ipoin];
		}
	}
	
	//switch order of element connectivity if(nnode==8)
	/* If elements are eight-node isoparametric elements,
	   rearrange the connecitivity list as required by vtk file format. */
	if(nnode == 8){
		for(int ielem=0;ielem<nelem;ielem++){
			for(int inode=0;inode<nnode;inode++){
				dummy[inode] = lnods[ielem][inode];
			}
			lnods[ielem][0] = dummy[0];
			lnods[ielem][1] = dummy[2];
			lnods[ielem][2] = dummy[4];
			lnods[ielem][3] = dummy[6];
			lnods[ielem][4] = dummy[1];
			lnods[ielem][5] = dummy[3];
			lnods[ielem][6] = dummy[5];
			lnods[ielem][7] = dummy[7];
		}
	}//end if
	
	//output to vtk file
	/* Write results in vtk format to be viewed by using Paraview.
	   If the deformed mesh configuration is desired, replace the coord array with
	   the dis_cord array which is calculated in lines "Add displacement component values,...." */
	istep = 1;
	write_vtk_fem(npoin,nelem,nnode,lnods,coord,istep,node_stres);
	
	return;
}