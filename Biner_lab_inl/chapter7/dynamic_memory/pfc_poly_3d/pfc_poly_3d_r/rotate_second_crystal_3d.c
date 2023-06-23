#include <stdio.h>
#include <stdlib.h> //rand()
#include <math.h> //M_PI

void rotate_second_crystal_3d(double theta_x, double theta_y, double theta_z,
	double con, int Nx1, int Ny1, int Nz1, double *den0, double *den1){
	
	double tmatx[3][3];
	double tmaty[3][3];
	double tmatz[3][3];
	
	double vect[3];
	//double vectx[3];
	double vecty[3];
	double vectz[3];
	double vect2[3];
	
	int ix;
	int iy;
	int iz;
	
	double pix=4.0*atan(1.0);
	
	//rotation matrix data (convert theta into radian)
	theta_x=theta_x*(pix/180.0);
	theta_y=theta_y*(pix/180.0);
	theta_z=theta_z*(pix/180.0);
	//
	//[Rx(theta_x)]^t (transposed matrix)
	tmatx[0][0] =  1.0;
	tmatx[0][1] =  0.0;
	tmatx[0][2] =  0.0;
	tmatx[1][0] =  0.0;
	tmatx[1][1] =  cos(theta_x);
	tmatx[1][2] =  sin(theta_x);
	tmatx[2][0] =  0.0;
	tmatx[2][1] = -sin(theta_x);
	tmatx[2][2] =  cos(theta_x);
	//
	//[Ry(theta_y)]^t (transposed matrix)
	tmaty[0][0] =  cos(theta_y);
	tmaty[0][1] =  0.0;
	tmaty[0][2] = -sin(theta_y);
	tmaty[1][0] =  0.0;
	tmaty[1][1] =  1.0;
	tmaty[1][2] =  0.0;
	tmaty[2][0] =  sin(theta_y);
	tmaty[2][1] =  0.0;
	tmaty[2][2] =  cos(theta_y);
	//
	//[Rz(theta_z)]^t (transposed matrix)
	tmatz[0][0] =  cos(theta_z);
	tmatz[0][1] =  sin(theta_z);
	tmatz[0][2] =  0.0;
	tmatz[1][0] = -sin(theta_z);
	tmatz[1][1] =  cos(theta_z);
	tmatz[1][2] =  0.0;
	tmatz[2][0] =  0.0;
	tmatz[2][1] =  0.0;
	tmatz[2][2] =  1.0;
	
	//rotate second crystal
	for(int i=0;i<Nx1;i++){
		for(int j=0;j<Ny1;j++){
			for(int k=0;k<Nz1;k++){
				//den1[i][j][k]=con;
				den1[i*Ny1*Nz1+j*Nz1+k]=con;
				
				vect[0]=i;
				vect[1]=j;
				vect[2]=k;
				
				//Rz
				for(int ii=0;ii<3;ii++){
					vectz[ii]=0.0;
					for(int jj=0;jj<3;jj++){
						vectz[ii]=vectz[ii]+tmatz[ii][jj]*vect[jj];
					}
				}
				//Ry
				for(int ii=0;ii<3;ii++){
					vecty[ii]=0.0;
					for(int jj=0;jj<3;jj++){
						vecty[ii]=vecty[ii]+tmaty[ii][jj]*vectz[jj];
					}
				}
				//Rx
				for(int ii=0;ii<3;ii++){
					vect2[ii]=0.0;
					for(int jj=0;jj<3;jj++){
						vect2[ii]=vect2[ii]+tmatx[ii][jj]*vecty[jj];
					}
				}
				
				ix=(int)(fmod(vect2[0],Nx1));
				iy=(int)(fmod(vect2[1],Ny1));
				iz=(int)(fmod(vect2[2],Nz1));
				
				//den1[i][j][k]=den0[ix][iy][iz];
				den1[i*Ny1*Nz1+j*Nz1+k]=den0[ix*Ny1*Nz1+iy*Nz1+iz];
			}
		}
	}
	
	return;
}