/* Prepare PFC bi-crystal */

#include <stdio.h>
#include <stdlib.h>
#include <math.h> //mod() and -lm

#define bicrystal 1
#define width 16.0
#define theta_rotate_second_crystal 15.0

//bicrystal==1 case (circle type)
#define radius_size_of_second_crystal 120.0

//bicrystal==2 case (side type)
#define left_side 120.0

void write_vtk_grid_values_2D(int nx, int ny, 
	double dx, double dy,
	int istep, double *data1, double *data2);

int main(){
	
	FILE *in=fopen("final_conf.out","r");
	FILE *out=fopen("bi_2r.inp","w");
	
	//simulation cell parameters
	int Nx=64;
	int Ny=64;
	//
	int Ntimes_x=8;
	int Ntimes_y=8;
	//
	int Nx1=(Nx*Ntimes_x);
	int Ny1=(Ny*Ntimes_y);
	
	//double data2[Nx1][Ny1];
	double *data2 = (double *)malloc(sizeof(double)*( Nx1*Ny1 ));
	//double   den[Nx1][Ny1];
	double *den   = (double *)malloc(sizeof(double)*( Nx1*Ny1 ));
	//double  den0[Nx1][Ny1];
	double *den0  = (double *)malloc(sizeof(double)*( Nx1*Ny1 ));
	//double  den1[Nx1][Ny1];
	double *den1  = (double *)malloc(sizeof(double)*( Nx1*Ny1 ));
	
	//The value of pi
	double pix=4.0*atan(1.0);
	
	//The distance between two grid points in x,y,z-direction
	double dx=pix/4.0;
	double dy=pix/4.0;
	
	double con=-0.285;
	
	double radius;
	
	double xlength;
	
	double theta;
	double tmatx[2][2];
	double vect[2];
	double vect2[2];
	
	int mx;
	int my;
	
	int ix;
	int iy;
	
	int ii;
	int jj;
	
	int iii0;
	int iii1;
	int iii2;
	
	int istep=1;
	
	//fscanf(in,"%5d %5d",&mx,&my);
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			iii0=i*Ny1+j;
			fscanf(in,"%5d %5d %lf",&mx,&my,&den0[iii0]);
		}
	}
	
	//replacate 64x64 grid to 512x512 grid
	for(int k=0;k<Ntimes_x;k++){
		for(int i=0;i<Nx;i++){
			ii=k*Nx+i;
			for(int m=0;m<Ntimes_y;m++){
				for(int j=0;j<Ny;j++){
					jj=m*Ny+j;
					iii0=i*Ny1+j;
					iii1=ii*Ny1+jj;
					den0[iii1]=den0[iii0];
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
				iii2=i*Ny1+j;
				den[iii2]=den0[iii2];
				xlength=sqrt( (i-Nx1/2)*(i-Nx1/2) + (j-Ny1/2)*(j-Ny1/2) );
				if(xlength<=radius){
					den[iii2]=con;
				}
			}
		}
		
		//rotate second crystal
		//theta=15.0*pix/180.0;
		theta=(theta_rotate_second_crystal)*pix/180.0;
		
		tmatx[0][0] =  cos(theta);
		tmatx[0][1] =  sin(theta);
		tmatx[1][0] = -sin(theta);
		tmatx[1][1] =  cos(theta);
		
		for(int i=0;i<Nx1;i++){
			for(int j=0;j<Ny1;j++){
				iii2=i*Ny1+j;
				den1[iii2]=con;
				
				vect[0]=i;
				vect[1]=j;
				
				for(int ii=0;ii<2;ii++){
					vect2[ii]=0.0;
					for(int jj=0;jj<2;jj++){
						vect2[ii]=vect2[ii]+tmatx[ii][jj]*vect[jj];
					}
				}
				
				ix=(int)(fmod(vect2[0],Nx1));
				iy=(int)(fmod(vect2[1],Ny1));
				
				//if(ix<=-1){
				//	ix=ix+Nx1;
				//}
				//if(ix>=Nx1){
				//	ix=ix-Nx1;
				//}
				//
				//if(iy<=-1){
				//	iy=iy+Ny1;
				//}
				//if(iy>=Ny1){
				//	iy=iy-Ny1;
				//}
				
				den1[iii2]=den0[ix*Ny1+iy];
			}
		}
		
		//size of second crystal
		//radius=120.0;
		radius=radius_size_of_second_crystal;
		
		for(int i=0;i<Nx1;i++){
			for(int j=0;j<Ny1;j++){
				xlength=sqrt( (i-Nx1/2)*(i-Nx1/2) + (j-Ny1/2)*(j-Ny1/2) );
				if(xlength<=radius){
					iii2=i*Ny1+j;
					den[iii2]=den1[iii2];
				}
			}
		}
	}//if
	
	if(bicrystal==2){
		//liquid region
		for(int i=0;i<Nx1;i++){
			for(int j=0;j<Ny1;j++){
				iii2=i*Ny1+j;
				den[iii2]=den0[iii2];
				if( (i>=(int)left_side) && (i<=(int)(Nx1-left_side)) ){
					den[iii2]=con;
				}
			}
		}
		
		//second crsytal
		
		//theta=30.0*pix/180.0;
		theta=(theta_rotate_second_crystal)*pix/180.0;
		
		tmatx[0][0] =  cos(theta);
		tmatx[0][1] =  sin(theta);
		tmatx[1][0] = -sin(theta);
		tmatx[1][1] =  cos(theta);
		
		for(int i=0;i<Nx1;i++){
			for(int j=0;j<Ny1;j++){
				iii2=i*Ny1+j;
				den1[iii2]=con;
				
				vect[0]=i;
				vect[1]=j;
				
				for(int ii=0;ii<2;ii++){
					vect2[ii]=0.0;
					for(int jj=0;jj<2;jj++){
						vect2[ii]=vect2[ii]+tmatx[ii][jj]*vect[jj];
					}
				}
				
				ix=(int)(fmod(vect2[0],Nx1));
				iy=(int)(fmod(vect2[1],Ny1));
				
				//if(ix<=-1){
				//	ix=ix+Nx1;
				//}
				//if(ix>=Nx1){
				//	ix=ix-Nx1;
				//}
				//
				//if(iy<=-1){
				//	iy=iy+Ny1;
				//}
				//if(iy>=Ny1){
				//	iy=iy-Ny1;
				//}
				
				den1[iii2]=den0[ix*Ny1+iy];
			}
		}
		
		//size of the second grain
		for(int i=(int)(left_side+width);i<(int)((Nx1-left_side)-width);i++){
			for(int j=0;j<Ny1;j++){
				iii2=i*Ny1+j;
				den[iii2]=den1[iii2];
			}
		}
	}//if
	
	//print output
	for(int i=0;i<Nx1;i++){
		for(int j=0;j<Ny1;j++){
			iii2=i*Ny1+j;
			fprintf(out,"%5d %5d %14.6e \n",i,j,den[iii2]);
		}
	}
	
	//graphic output
	//istep=1;
	for(int i=0;i<Nx1;i++){
		for(int j=0;j<Ny1;j++){
			iii2=i*Ny1+j;
			data2[iii2]=0.0;
		}
	}
	
	write_vtk_grid_values_2D(Nx1,Ny1,dx,dy,istep,den,data2);
	
	fclose(in);
	fclose(out);
}
