/* This function evaluates the values of shape functions and
   their derivatives in local coordinates. */

/* Variable and array list
  exisp: x-coordinate of sampling point in local coordinates (double)
  etasp: y-coordinate of sampling point in local coordinates (double)
  nnode: Number of nodes per element (int)
  shape[nnode]: Shape functions values for each node in an element (double)
  deriv[ndime][nnode]: Derivatives of shape functions for
    each local coordinate direction. (double)
*/

void sfr2_2d(double exisp, double etasp, int nnode,
	double *shape, double *deriv){
	
	double s, t, p;
	
	if(nnode == 3){
		//3 nodes elements
		/* The values of shape functions and their derivatives are
		   calculated for three-node isoparametric elements. */
		
		s = exisp;
		t = etasp;
		p = 1.0 - s - t;
		
		shape[0] = p;
		shape[1] = s;
		shape[2] = t;
		
		deriv[0*nnode+0] = -1.0;
		deriv[0*nnode+1] =  1.0;
		deriv[0*nnode+2] =  0.0;
		
		deriv[1*nnode+0] = -1.0;
		deriv[1*nnode+1] =  0.0;
		deriv[1*nnode+2] =  1.0;
	}
	
	double st;
	
	if(nnode == 4){
		//4 nodes elements
		/* The values of shape functions and their derivatives are
		   calculated for four-node isoparametric elements. */
		
		s = exisp;
		t = etasp;
		st = s * t;
		
		shape[0] = (1.0 - t - s + st)*0.25;
		shape[1] = (1.0 - t + s - st)*0.25;
		shape[2] = (1.0 + t + s + st)*0.25;
		shape[3] = (1.0 + t - s - st)*0.25;
		
		deriv[0*nnode+0] = (-1.0 + t)*0.25;
		deriv[0*nnode+1] = ( 1.0 - t)*0.25;
		deriv[0*nnode+2] = ( 1.0 + t)*0.25;
		deriv[0*nnode+3] = (-1.0 - t)*0.25;
		
		deriv[1*nnode+0] = (-1.0 + s)*0.25;
		deriv[1*nnode+1] = (-1.0 - s)*0.25;
		deriv[1*nnode+2] = ( 1.0 + s)*0.25;
		deriv[1*nnode+3] = ( 1.0 - s)*0.25;
	}
	
	double s2, t2;
	double ss, tt;
	double sst, stt;
	double st2;
	
	if(nnode == 8){
		//8 nodes elements
		/* The values of shape functions and their derivatives are
		   calculated for eight-node isoparametric elements. */
		
		s = exisp;
		t = etasp;
		s2 = 2.0 * s;
		t2 = 2.0 * t;
		ss = s * s;
		tt = t * t;
		st = s * t;
		sst = s * s * t;
		stt = s * t * t;
		st2 = 2.0 * s * t;
		
		shape[0] = (-1.0 + st + ss + tt - sst - stt)*0.25;
		shape[1] = ( 1.0 -  t - ss + tt + sst      )*0.25;
		shape[2] = (-1.0 - st + ss + tt - sst + stt)*0.25;
		shape[3] = ( 1.0 +  s      - tt       - stt)*0.25;
		shape[4] = (-1.0 + st + ss + tt + sst + stt)*0.25;
		shape[5] = ( 1.0 +  t - ss      - sst - stt)*0.25;
		shape[6] = (-1.0 - st + ss + tt + sst - stt)*0.25;
		shape[7] = ( 1.0 -  s +    - tt       + stt)*0.25;
		
		deriv[0*nnode+0] = ( t + s2 - st2 - tt)*0.25;
		deriv[0*nnode+1] = (-s + st);
		deriv[0*nnode+2] = (-t + s2 - st2 + tt)*0.25;
		deriv[0*nnode+3] = ( 1.0 - tt)*0.5;
		deriv[0*nnode+4] = ( t + s2 + st2 + tt)*0.25;
		deriv[0*nnode+5] = (-s - st);
		deriv[0*nnode+6] = (-t + s2 + st2 -tt)*0.25;
		deriv[0*nnode+7] = (-1.0 + tt)*0.5;
		
		deriv[1*nnode+0] = ( s + t2 - ss  - st2)*0.25;
		deriv[1*nnode+1] = (-1.0 + ss)*0.5;
		deriv[1*nnode+2] = (-s + t2 - ss  + st2)*0.25;
		deriv[1*nnode+3] = (-t - st);
		deriv[1*nnode+4] = ( s + t2 + ss  + st2)*0.25;
		deriv[1*nnode+5] = ( 1.0 - ss)*0.5;
		deriv[1*nnode+6] = (-s + t2 + ss  - st2)*0.25;
		deriv[1*nnode+7] = (-t + st);
	}
	
	return;
}