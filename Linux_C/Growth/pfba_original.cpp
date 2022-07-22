#define MESH 300
#define PAT 3.14159265
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
//#include <conio.h>

void init_cond(float temp[][MESH], float phase[][MESH]);
void phase_change(float temp[][MESH], float phase[][MESH], 
	              float tau, float eps, 
	              float alpha, float gamma, float del_t,
                  float del_x, float latent_heat);
void write_to_file(float temp[][MESH], float phase[][MESH], int step);

FILE *fp;

int main(void)
{
	float temp[MESH][MESH];
	float phase[MESH][MESH];
	float tau = 0.0003;      
	float eps = 0.01;        
	float alpha = 0.9;       
	float gamma = 10.0;      
	float latent_heat = 1.0; 
	float del_t = 0.0002;    
	float t_max = 9.6;       
	float del_x = 0.03;      
	int step, max_step, w_times = 10;
	
	if((fp = fopen("phase4.dat","w")) == NULL){
		printf("file open error!");
		fclose(fp);
		//return;
	}
	
	max_step = (int)(t_max / del_t);
	
	init_cond(temp, phase);
	write_to_file(temp, phase, 0);
	
	for( step = 1; step <= max_step; step++ ){
		printf("step = %1d \n", step);
		phase_change(temp, phase, tau, eps, alpha, gamma, del_t, del_x, latent_heat);
		if( ( step % (max_step / w_times)) == 0 ){
			write_to_file(temp, phase, step);
		}
	}
	fclose(fp);
}

void init_cond(float temp[][MESH], float phase[][MESH])
{
	int i, j;
	for(i = 0; i < MESH; i++){
		for(j = 0; j < MESH; j++){
			if( i == 0 | i == MESH - 1 | j == MESH - 1){
				temp[i][j] = 0.0; phase[i][j] = 1.0;
				// wall is cool and solid state
			} else {
				temp[i][j] = 1.0; phase[i][j] = 0.0;
				// other sites are in T_equilibrium and liquid
			}
		}
	}
}

void phase_change(float temp[][MESH], float phase[][MESH], 
	              float tau, float eps, 
	              float alpha, float gamma, float del_t,
	              float del_x, float latent_heat)
{
	int i, j, i_left, i_right, j_up, j_down;
	float const_e, const_p, const_t, m_t, d_t;
	float p_left, p_right, p_up, p_down, p_ij;
	float phase_change[MESH][MESH];
	
	const_e = eps*eps * del_t / tau / del_x / del_x;
	const_p = del_t / tau;
	d_t = del_t / del_x / del_x;
	const_t = alpha / PAT;
	
	for(i = 0; i < MESH; i++){
		for(j = 0; j < MESH; j++){
			p_ij = phase[i][j];
			m_t = const_t * (float) atan( (double)(gamma*(1.0 - temp[i][j])) );
			i_left = i - 1; i_right = i + 1;
			j_up = j + 1; j_down = j - 1;
			
			if( i_left == -1 ){
				phase[i][j] = 1.0; continue;
			} else {
				p_left = phase[i_left][j];
			}
			
			if( i_right == MESH ){
				phase[i][j] = 1.0; continue;
			} else {
				p_right = phase[i_right][j];
			}
			
			if( j_up == MESH ){
				phase[i][j] = 1.0; continue;
			} else {
				p_up = phase[i][j_up];
			}
			
			if( j_down == -1 ){
				phase[i][j] = 1.0; continue;
			} else {
				p_down = phase[i][j_down];
			}
			
			phase_change[i][j] = const_e * ( p_right + p_left + p_up + p_down - 4.0*p_ij) 
			                   + const_p*p_ij*(1.0 - p_ij)*(p_ij - 0.5 + m_t);
			phase[i][j] = phase_change[i][j] + p_ij;
		}
	}
	
	for( i = 0; i < MESH; i++){
		for( j = 0; j < MESH; j++){
			p_ij = temp[i][j];
			i_left = i - 1; i_right = i + 1;
			j_up = j + 1; j_down = j - 1;
			
			if( i_left == -1 ){
				temp[i][j] = 0; continue;
			} else {
				p_left = temp[i_left][j];
			}
			
			if( i_right == MESH ){
				temp[i][j] = 0; continue;
			} else {
				p_right = temp[i_right][j];
			}
			
			if( j_up == MESH ){
				temp[i][j] = 0; continue;
			} else {
				p_up = temp[i][j_up];
			}
			
			if( j_down == -1){
				temp[i][j] = 0; continue;
			} else {
				p_down = temp[i][j_down];
			}
			
			temp[i][j] = temp[i][j] + d_t*( p_right + p_left + p_up + p_down - 4.0*temp[i][j] )
			           + latent_heat*phase_change[i][j];
		}
	}
}

void write_to_file(float temp[][MESH], float phase[][MESH], int step)
{
	int i, j;
	
	fprintf( fp, "%1d \n", step );
	for( i = 0; i < MESH; i++ ){
		for( j = 0; j < MESH; j++){
			fprintf( fp, " %f", temp[i][j] );
		}
		fprintf( fp, "\n" );
	}
	
	fprintf( fp, "\n" );
	
	for( i = 0; i < MESH; i++ ){
		for( j = 0; j < MESH; j++ ){
			fprintf( fp, " %f", phase[i][j]);
		}
		fprintf( fp, "\n" );
	}
	fprintf( fp, "\n" );
	fprintf( fp, "\n" );
}