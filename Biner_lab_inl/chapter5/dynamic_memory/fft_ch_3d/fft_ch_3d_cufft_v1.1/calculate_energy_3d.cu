/* This function calculates the total bulk energy of
   the system. The results are output to file
   time_energy.out to generate time vs. energy plots */

/* Variable and array list
  Nx: Number of grid points in the x-direction
  Ny: Number of grid points in the y-direction
  Nz: Number of grid points in the z-direction
  grad_coef: Gradient energy coefficient
  energy: total bulk energy value
  con[Nx][Ny][Nz]: Concentration filed
*/

float calculate_energy_3d(int Nx, int Ny, int Nz, float *con, float grad_coef){
	
	//initialize
	float energy=0.0;
	
	int ijk;
	int ipjk;
	int ijpk;
	int ijkp;
	
	int ip,jp,kp;
	
	//calculate Eq.4.19
	for(int i=0;i<(Nx-1);i++){
		ip = i + 1;
		for(int j=0;j<(Ny-1);j++){
			jp = j + 1;
			for(int k=0;k<(Nz-1);k++){
				kp = k + 1;
				//
				ijk=(i*Ny+j)*Nz+k;
				ipjk=(ip*Ny+j)*Nz+k;
				ijpk=(i*Ny+jp)*Nz+k;
				ijkp=(i*Ny+j)*Nz+kp;
				//
				energy = energy
					   +con[ijk]*con[ijk]*(1.0-con[ijk])*(1.0-con[ijk])
					   +0.5*grad_coef*((con[ipjk]-con[ijk])*(con[ipjk]-con[ijk])
					   				  +(con[ijpk]-con[ijk])*(con[ijpk]-con[ijk])
					   				  +(con[ijkp]-con[ijk])*(con[ijkp]-con[ijk])
					   				  );
				//
			}//end for(k
		}//end for(j
	}//end for(i
	
	return energy;
}