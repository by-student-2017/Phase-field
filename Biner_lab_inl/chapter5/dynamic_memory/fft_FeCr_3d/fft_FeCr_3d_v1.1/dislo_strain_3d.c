/* This function introduces the eigenstrain distribution of
   dislocations to the simulation cell. idislo=1 is for
   a dipole and idislo=2 is for dislocation array */

/* Variable and array list
  Nx: Number of grid points in the x-direction
  Ny: Number of grid points in the y-direction
  Nz: Number of grid points in the z-direction
  idislo: idislo=1 is for dislocsation dipole, and 
          idislo=2 is for dislocation array
  ed11[Nx][Ny][Nz]: Eigenstrain component of dislocations (see Eq.5.54)
  ed22[Nx][Ny][Nz]: Eigenstrain component of dislocations
  ed12[Nx][Ny][Nz]: Eigenstrain component of dislocations
*/

void dislo_strain_3d(int Nx, int Ny, int Nz, int idislo,
	double *ed11, double *ed22, double *ed33,
	double *ed12, double *ed23, double *ed13){
	
	int ijk;
	
	//Initialize eigenstrain field of dislocations
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			for(int k=0;k<Nz;k++){
				ijk=(i*Ny+j)*Nz+k;
				//----- ----- -----
				ed11[ijk] = 0.0;
				ed22[ijk] = 0.0;
				ed33[ijk] = 0.0;
				//
				ed12[ijk] = 0.0;
				ed23[ijk] = 0.0;
				ed13[ijk] = 0.0;
				//----- ----- -----
			}
		}
	}
	
	int Ny2=(int)Ny/2;
	int Nz2=(int)Nz/2;
	
	// Begin and end grid points of dislocation dipoles
	int ndipoles=(34-1);
	int ndipolee=(94-1);
	
	// Introduce eigenstrain values of dislocation dipole into the simulation cell
	if(idislo==1){
		for(int i=0;i<Nx;i++){
			//if(i>=34 && i<=94){
			if(i>=ndipoles && i<=ndipolee){
				ed12[(i*Ny+Ny2)*Nz+Nz2]=5.0e-3;
			}
		}
	}//end if
	
	int ndis=11;
	int jj=0;
	
	/* Introduce eigenstrain values of a dislocation array,
	   12 grid points apart in in the y-direction, into the simulation cell */
	if(idislo==2){
		//
		for(int idis=0;idis<ndis;idis++){
			jj=jj+12;
			for(int i=0;i<Nx;i++){
				for(int j=0;j<Ny;j++){
					//if(i>=34 && i<94){
					if(i>=ndipoles && i<ndipolee){
						if(j==jj){
							ijk=(i*Ny+j)*Nz+Nz2; //k=Nz2
							ed12[ijk]=5.0e-3;
						}
					}
				}//j
			}//i
		}//idis
		//
	}//end if
	
	return;
}