/* Prepare PFC bi-crystal */

#include <stdio.h>
#include <math.h> //mod() and -lm

#define bicrystal 1
#define width 16.0
#define theta_rotate_second_crystal 15.0

//bicrystal==1 case (circle type)
#define radius_size_of_second_crystal 120.0

//bicrystal==2 case (side type)
#define left_side 120.0

#define Nx 64
#define Ny 64

#define Ntimes_x 8
#define Ntimes_y 8

#define Nx1 (Nx*Ntimes_x) //e.g., 64*8=512
#define Ny1 (Ny*Ntimes_y) //e.g., 64*8=512

	double data2[Nx1][Ny1];
	double   den[Nx1][Ny1];
	double  den0[Nx1][Ny1];
	double  den1[Nx1][Ny1];

void write_vtk_grid_values_2D();

int main(){
	
	FILE *in=fopen("final_conf.out","r");
	FILE *out=fopen("bi_2r.inp","w");
	
	//simulation cell parameters
	//int Nx=64;
	//int Ny=64;
	
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
	
	int istep=1;
	
	//fscanf(in,"%5d %5d",&mx,&my);
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			fscanf(in,"%5d %5d %lf",&mx,&my,&den0[i][j]);
		}
	}
	
	//replacate 64x64 grid to 512x512 grid
	for(int k=0;k<Ntimes_x;k++){
		for(int i=0;i<Nx;i++){
			ii=k*Nx+i;
			for(int m=0;m<Ntimes_y;m++){
				for(int j=0;j<Ny;j++){
					jj=m*Ny+j;
					den0[ii][jj]=den0[i][j];
				}
			}
		}
	}
	//
	//for(int i=0;i<Nx;i++){
	//	for(int m=0;m<Ntimes_y;m++){
	//		for(int j=0;j<Ny;j++){
	///			jj=m*Ny+j;
	//			den0[i][jj]=den0[i][j];
	//		}
	//	}
	//}
	//
	//for(int k=0;k<Ntimes_x;k++){
	//	for(int i=0;i<Nx;i++){
	//		ii=Nx*k+i;
	//		for(int j=0;j<Ny;j++){
	//			den0[ii][j]=den0[i][j];
	//		}
	//	}
	//}
	
	//generate bicrystal
	//for(int i=0;i<Nx1;i++){
	//	for(int j=0;j<Ny1;j++){
	//		den[i][j]=con;
	//	}
	//}
	
	if(bicrystal==1){
		//liquid region
		//radius=142.0;
		radius=(radius_size_of_second_crystal+width);
		
		for(int i=0;i<Nx1;i++){
			for(int j=0;j<Ny1;j++){
				den[i][j]=den0[i][j];
				xlength=sqrt( (i-Nx1/2)*(i-Nx1/2) + (j-Ny1/2)*(j-Ny1/2) );
				if(xlength<=radius){
					den[i][j]=con;
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
				den1[i][j]=con;
				
				vect[0]=i;
				vect[1]=j;
				
				for(int ii=0;ii<2;ii++){
					vect2[ii]=0.0;
					for(int jj=0;jj<2;jj++){
						vect2[ii]=vect2[ii]+tmatx[ii][jj]*vect[jj];
					}
				}
				
				ix=(int)vect2[0];
				iy=(int)vect2[1];
				
				if(ix<=-1){
					ix=ix+Nx1;
				}
				if(ix>=Nx1){
					ix=ix-Nx1;
				}
				//
				if(iy<=-1){
					iy=iy+Ny1;
				}
				if(iy>=Ny1){
					iy=iy-Ny1;
				}
				
				den1[i][j]=den0[ix][iy];
			}
		}
		
		//size of second crystal
		//radius=120.0;
		radius=radius_size_of_second_crystal;
		
		for(int i=0;i<Nx1;i++){
			for(int j=0;j<Ny1;j++){
				xlength=sqrt( (i-Nx1/2)*(i-Nx1/2) + (j-Ny1/2)*(j-Ny1/2) );
				if(xlength<=radius){
					den[i][j]=den1[i][j];
				}
			}
		}
	}//if
	
	if(bicrystal==2){
		//liquid region
		for(int i=0;i<Nx1;i++){
			for(int j=0;j<Ny1;j++){
				den[i][j]=den0[i][j];
				if( (i>=(int)left_side) && (i<=(int)(Nx1-left_side)) ){
					den[i][j]=con;
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
				den1[i][j]=con;
				
				vect[0]=i;
				vect[1]=j;
				
				for(int ii=0;ii<2;ii++){
					vect2[ii]=0.0;
					for(int jj=0;jj<2;jj++){
						vect2[ii]=vect2[ii]+tmatx[ii][jj]*vect[jj];
					}
				}
				
				ix=(int)vect2[0];
				iy=(int)vect2[1];
				
				if(ix<=-1){
					ix=ix+Nx1;
				}
				if(ix>=Nx1){
					ix=ix-Nx1;
				}
				//
				if(iy<=-1){
					iy=iy+Ny1;
				}
				if(iy>=Ny1){
					iy=iy-Ny1;
				}
				
				den1[i][j]=den0[ix][iy];
			}
		}
		
		//size of the second grain
		for(int i=(int)(left_side+width);i<(int)((Nx1-left_side)-width);i++){
			for(int j=0;j<Ny1;j++){
				den[i][j]=den1[i][j];
			}
		}
	}//if
	
	//print output
	for(int i=0;i<Nx1;i++){
		for(int j=0;j<Ny1;j++){
			fprintf(out,"%5d %5d %14.6e \n",i,j,den[i][j]);
		}
	}
	
	//graphic output
	//istep=1;
	for(int i=0;i<Nx1;i++){
		for(int j=0;j<Ny1;j++){
			data2[i][j]=0.0;
		}
	}
	
	write_vtk_grid_values_2D(Nx1,Ny1,dx,dy,istep,den,data2);
	
	fclose(in);
	fclose(out);
}
