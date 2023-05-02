/* This function introduces the eigenstrain distribution of
   dislocations to the simulation cell. idislo=1 is for
   a dipole and idislo=2 is for dislocation array */

#include <stdlib.h> //rand()
#include <fftw3.h>

/* Variable and array list
  Nx: Number of grid points in the x-direction
  Ny: Number of grid points in the y-direction
  idislo: idislo=1 is for dislocsation dipole, and 
          idislo=2 is for dislocation array
  ed11[Nx][Ny]: Eigenstrain component of dislocations (see Eq.5.54)
  ed22[Nx][Ny]: Eigenstrain component of dislocations
  ed12[Nx][Ny]: Eigenstrain component of dislocations
*/

void dislo_strain_2d(int Nx, int Ny, int idislo,
	fftw_complex *ed11, fftw_complex *ed22, fftw_complex *ed12){
	
	int ii;
	
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			ii=i*Ny+j;
			//
			ed11[ii][0] = 0.0;
			ed11[ii][1] = 0.0;
			//
			ed22[ii][0] = 0.0;
			ed22[ii][1] = 0.0;
			//
			ed12[ii][0] = 0.0;
			ed12[ii][1] = 0.0;
		}
	}
	
	int Ny2=(int)Ny/2;
	
	int ndipoles=(34-1);
	int ndipolee=(94-1);
	
	if(idislo==1){
		for(int i=0;i<Nx;i++){
			//if(i>=34 && i<=94){
			if(i>=ndipoles && i<=ndipolee){
				ed12[i*Ny+Ny2][0]=5.0e-3;
				ed12[i*Ny+Ny2][1]=0.0;
			}
		}
	}//end if
	
	int ndis=12;
	int jj=0;
	
	if(idislo==2){
		//
		for(int idis=0;idis<ndis;idis++){
			//jj=jj+12;
			jj=jj+ndis;
			for(int i=0;i<Nx;i++){
				for(int j=0;j<Ny;j++){
					//if(i>=34 && i<94){
					if(i>=ndipoles && i<ndipolee){
						if(j==jj){
							ii=i*Ny+j;
							ed12[ii][0]=5.0e-3;
							ed12[ii][1]=0.0;
						}
					}
				}//j
			}//i
		}//idis
		//
	}//end if
	
	return;
}