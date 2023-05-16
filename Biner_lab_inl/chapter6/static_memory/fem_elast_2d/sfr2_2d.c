/* This function evaluates the values of shape functions and
   their derivatives in local coordinates. */

/* Variable and array list
  exisp: x-coordinate of sampling point in local coordinates
  etasp: y-coordinate of sampling point in local coordinates
  nnode: Number of nodes per element
  shape[nnode]: Shape functions values for each node in an element
  deriv[ndime][nnode]: Derivatives of shape functions for
    each local coordinate direction.
*/

void sfr2_2d(double exisp, double etasp, int nnode,
	double *shape, double *deriv){
	
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
		
		deriv[0][0] = -1.0;
		deriv[0][1] =  1.0;
		deriv[0][2] =  0.0;
		
		deriv[1][0] = -1.0;
		deriv[1][1] =  0.0;
		deriv[1][2] =  1.0;
	}
	
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
		
		deriv[0][0] = (-1.0 + t)*0.25;
		deriv[0][1] = ( 1.0 - t)*0.25;
		deriv[0][2] = ( 1.0 + t)*0.25;
		deriv[0][3] = (-1.0 - t)*0.25;
		
		deriv[1][0] = (-1.0 + s)*0.25;
		deriv[1][1] = (-1.0 - s)*0.25;
		deriv[1][2] = ( 1.0 + s)*0.25;
		deriv[1][3] = ( 1.0 - s)*0.25;
	}
	
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
		
		deriv[0][0] = ( t + s2 - st2 - tt)*0.25;
		deriv[0][1] = (-s + st);
		deriv[0][2] = (-t + s2 - st2 + tt)*0.25;
		deriv[0][3] = ( 1.0 - tt)*0.5;
		deriv[0][4] = ( t + s2 + st2 + tt)*0.25;
		deriv[0][5] = (-s - st);
		deriv[0][6] = (-t + s2 + st2 -tt)*0.25;
		deriv[0][7] = (-1.0 + tt)*0.5;
		
		deriv[1][0] = ( s + t2 - ss2 - st2)*0.25;
		deriv[1][1] = (-1.0 + ss)*0.5;
		deriv[1][2] = (-s + t2 - ss  + st2)*0.25;
		deriv[1][3] = (-t - st);
		deriv[1][4] = ( s + t2 + ss  + st2)*0.25;
		deriv[1][5] = ( 1.0 - ss)*0.5;
		deriv[1][6] = (-s + t2 + ss  - st2)*0.25;
		deriv[1][7] = (-t + st);
	}
	
	return;
}