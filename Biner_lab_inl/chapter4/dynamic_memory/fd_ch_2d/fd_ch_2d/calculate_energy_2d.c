/* This function calculates the total bulk energy of
   the system. The results are output to file
   time_energy.out to generate time vs. energy plots */

/* Variable and array list
  Nx: Number of grid points in the x-direction
  Ny: Number of grid points in the y-direction
  grad_coef: Gradient energy coefficient
  energy: total bulk energy value
  con[Nx][Ny]: Concentration filed
*/

double calculate_energy_2d(int Nx, int Ny, double *con, double grad_coef){
	
	//initialize
	double energy=0.0;
	
	int ij;
	int ipj;
	int ijp;
	
	int ip,jp;
	
	//calculate Eq.4.19
	for(int i=0;i<(Nx-1);i++){
		ip = i + 1;
		for(int j=0;j<(Ny-1);j++){
			jp = j + 1;
			//
			ij=i*Ny+j;
			ipj=ip*Ny+j;
			ijp=i*Ny+jp;
			//
			energy = energy
				   +con[ij]*con[ij]*(1.0-con[ij])*(1.0-con[ij])
				   +0.5*grad_coef*((con[ipj]-con[ij])*(con[ipj]-con[ij])
				   				  +(con[ijp]-con[ij])*(con[ijp]-con[ij])
				   				  );
		}
	}
	
	return energy;
}