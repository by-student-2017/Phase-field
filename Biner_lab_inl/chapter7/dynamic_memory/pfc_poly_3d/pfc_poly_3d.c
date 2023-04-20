/* Prepare PFC bi-crystal */

#include <stdio.h>
#include <stdlib.h>
#include <math.h> //mod() and -lm

#define bicrystal 2
#define width 10.0

//[r']=[Rx(theta_x)][Ry(theta_y)][Rz(theta_z)][r]
#define theta_x  0.0 //rotate_second_crystal
#define theta_y  0.0 //rotate_second_crystal
#define theta_z 30.0 //rotate_second_crystal

//bicrystal==1 case (circle type)
#define radius_size_of_second_crystal 24.0

//bicrystal==2 case (side type)
#define left_side 20.0

void rotate_second_crystal_3d(double theta_xx, double theta_yy, double theta_zz,
	double con, int Nx1, int Ny1, int Nz1, double *den0, double *den1);

void write_vtk_grid_values_3D(int nx, int ny, int nz, 
	double dx, double dy, double dz,
	int istep, double *data1, double *data2);

int main(){
	
	FILE *in=fopen("final_conf.out","r");
	FILE *out=fopen("bi_2r_3d.inp","w");
	
	//simulation cell parameters
	int Nx=32;
	int Ny=32;
	int Nz=32;
	//
	int Ntimes_x=3;
	int Ntimes_y=3;
	int Ntimes_z=3;
	//
	int Nx1=(Nx*Ntimes_x);
	int Ny1=(Ny*Ntimes_y);
	int Nz1=(Nz*Ntimes_z);
	
	//double  den0[Nx1][Ny1][Nz1];
	double *den0  = (double *)malloc(sizeof(double)*( Nx1*Ny1*Nz1 ));
	//double  den1[Nx1][Ny1][Nz1];
	double *den1  = (double *)malloc(sizeof(double)*( Nx1*Ny1*Nz1 ));
	//double   den[Nx1][Ny1][Nz1];
	double *den   = (double *)malloc(sizeof(double)*( Nx1*Ny1*Nz1 ));
	//double data2[Nx1][Ny1][Nz1];
	double *data2 = (double *)malloc(sizeof(double)*( Nx1*Ny1*Nz1 ));
	
	//The value of pi
	double pix=4.0*atan(1.0);
	
	//The distance between two grid points in x,y,z-direction
	double dx=pix/4.0;
	double dy=pix/4.0;
	double dz=pix/4.0;
	
	double con=-0.285;
	
	double radius;
	
	double xlength;
	
	int mx;
	int my;
	int mz;
	
	int ii;
	int jj;
	int kk;
	
	int iii0;
	int iii1;
	int iii2;
	
	int istep=1;
	
	//fscanf(in,"%5d %5d",&mx,&my);
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			for(int k=0;k<Nz;k++){
				iii0=i*Ny1*Nz1+j*Nz1+k;
				fscanf(in,"%5d %5d %5d %lf",&mx,&my,&mz,&den0[iii0]);
			}
		}
	}
	
	//replacate 32x32x32 grid to 96x96x96 grid
	for(int l=0;l<Ntimes_x;l++){
		for(int i=0;i<Nx;i++){
			ii=l*Nx+i;
			for(int m=0;m<Ntimes_y;m++){
				for(int j=0;j<Ny;j++){
					jj=m*Ny+j;
					for(int n=0;n<Ntimes_z;n++){
						for(int k=0;k<Nz;k++){
							kk=n*Nz+k;
							iii0=i*Ny1*Nz1+j*Nz1+k;
							iii1=ii*Ny1*Nz1+jj*Nz1+kk;
							den0[iii1]=den0[iii0];
						}
					}
				}
			}
		}
	}
	
	if(bicrystal==1){
		//liquid region
		//radius=142.0;
		radius=(radius_size_of_second_crystal+width);
		
		for(int i=0;i<Nx1;i++){
			for(int j=0;j<Ny1;j++){
				for(int k=0;k<Nz1;k++){
					iii2=i*Ny1*Nz1+j*Nz1+k;
					den[iii2]=den0[iii2];
					xlength=sqrt( (i-Nx1/2)*(i-Nx1/2)
								+ (j-Ny1/2)*(j-Ny1/2)
								+ (k-Nz1/2)*(k-Nz1/2) );
					if(xlength<=radius){
						den[iii2]=con;
					}
				}
			}
		}
		
		//rotate second crystal
		rotate_second_crystal_3d(theta_x, theta_y, theta_z,
			con, Nx1, Ny1, Nz1, den0, den1);
		
		//size of second crystal
		//radius=120.0;
		radius=radius_size_of_second_crystal;
		
		for(int i=0;i<Nx1;i++){
			for(int j=0;j<Ny1;j++){
				for(int k=0;k<Nz1;k++){
					xlength=sqrt( (i-Nx1/2)*(i-Nx1/2)
								+ (j-Ny1/2)*(j-Ny1/2)
								+ (k-Nz1/2)*(k-Nz1/2) );
					if(xlength<=radius){
						iii2=i*Ny1*Nz1+j*Nz1+k;
						den[iii2]=den1[iii2];
					}
				}
			}
		}
	}//if
	
	if(bicrystal==2){
		//liquid region
		for(int i=0;i<Nx1;i++){
			for(int j=0;j<Ny1;j++){
				for(int k=0;k<Nz1;k++){
					iii2=i*Ny1*Nz1+j*Nz1+k;
					den[iii2]=den0[iii2];
					if( (i>=(int)left_side) && (i<=(int)(Nx1-left_side)) ){
						den[iii2]=con;
					}
				}
			}
		}
		
		//rotate second crystal
		rotate_second_crystal_3d(theta_x, theta_y, theta_z,
			con, Nx1, Ny1, Nz1, den0, den1);
		
		//size of the second grain
		for(int i=(int)(left_side+width);i<(int)((Nx1-left_side)-width);i++){
			for(int j=0;j<Ny1;j++){
				for(int k=0;k<Nz1;k++){
					iii2=i*Ny1*Nz1+j*Nz1+k;
					den[iii2]=den1[iii2];
				}
			}
		}
	}//if
	
	//print output
	for(int i=0;i<Nx1;i++){
		for(int j=0;j<Ny1;j++){
			for(int k=0;k<Nz1;k++){
				iii2=i*Ny1*Nz1+j*Nz1+k;
				fprintf(out,"%5d %5d %5d %14.6e \n",i,j,k,den[iii2]);
			}
		}
	}
	
	//graphic output
	//istep=1;
	for(int i=0;i<Nx1;i++){
		for(int j=0;j<Ny1;j++){
			for(int k=0;k<Nz1;k++){
				iii2=i*Ny1*Nz1+j*Nz1+k;
				data2[iii2]=0.0;
			}
		}
	}
	
	write_vtk_grid_values_3D(Nx1,Ny1,Nz1,dx,dy,dz,istep,den,data2);
	
	fclose(in);
	fclose(out);
}
