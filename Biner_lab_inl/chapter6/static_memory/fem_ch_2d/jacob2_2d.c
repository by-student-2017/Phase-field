/* This function evaluates the Cartesian shape function derivatives,
   determinant of the Jacobian for the numerical integration of
   area/volume of elements, and the Cartesian coordinates of
   the integration points within elements. */

/* Variable and array list
  ielem: Current element number (int)
  kgasp: The current integration point number (int)
  ndime: Number of global coordinate components (int)
  nnode: Number of nodes per element (int)
  djacb: The determinate of the Jacobian matrix sampled within
    the element (ielem) at the integration points and
    used in calculation of area/volume of the elements (double)
  shape[nnode]: Shape function values (double)
  elcod[ndime][nnode]: Global coordinates of current element (ielem) nodes. (double)
  deriv[ndime][nnode]: Derivatives of shape functions is local coordinates. (double)
  gpcod[ndime][ngaus*ngaus]: Cartesian coordinates of the integration points within the elemet (ielem). (double)
*/

void jacob2_2d(int ielem, 
	double *elcod,
	int kgasp,
	double *shape, double *deriv,
	int nnode, int ndime,
	double *cartd, double djacb, double *gpcod,
	int ngaus){
	
	//Gauss point coordinates
	/* Calculate the Cartesian coordinates of
	   the integration points within the elements. */
	for(int idime=0;idime<ndime;idime++){
		gpcod[idime*(ngaus*ngaus)+kgasp]=0.0;
		for(int inode=0;inode<nnode;inode++){
			gpcod[idime*(ngaus*ngaus)+kgasp] = gpcod[idime*(ngaus*ngaus)+kgasp]
								+ elcod[idime*nnode+inode] * shape[inode];
		}
	}
	
	//Jacobian
	double xjacm[ndime][ndime];
	// Form the Jacobian Matrix (Eq.6.5)
	for(int idime=0;idime<ndime;idime++){
		for(int jdime=0;jdime<ndime;jdime++){
			xjacm[idime][jdime]=0.0;
			for(int inode=0;inode<nnode;inode++){
				xjacm[idime][jdime] = xjacm[idime][jdime]
									+ deriv[idime*nnode+inode] * elcod[jdime*nnode+inode];
			}
		}
	}
	
	/* Calculate the determinant of the Jacobian materix.
	   If its values is zero or negative, terminate the overall solution. */
	djacb = xjacm[0][0] * xjacm[1][1] - xjacm[0][1] * xjacm[1][0];
	
	if(djacb <= 0.0){
		printf("Program Terminated \n");
		printf("Zero or negative area for element: %5d \n",ielem);
		printf("Program terminated zero or negative area \n");
	}
	
	//Cartesian derivatives
	// Calculate the inverse of the Jacobian matrix (Eq.6.6)
	double xjaci[ndime][ndime];
	xjaci[0][0] = xjacm[1][1]/djacb;
	xjaci[1][1] = xjacm[0][0]/djacb;
	xjaci[0][1] = xjacm[0][1]/djacb;
	xjaci[1][0] = xjacm[1][0]/djacb;
	
	/* Calculate Cartesian derivatives of
	   shape functions using chain rule (Eq.6.3). */
	for(int idime=0;idime<ndime;idime++){
		for(int inode=0;inode<nnode;inode++){
			cartd[idime*nnode+inode]=0.0;
			for(int jdime=0;jdime<ndime;jdime++){
				cartd[idime*nnode+inode] = cartd[idime*nnode+inode]
										 + xjaci[idime][jdime] * deriv[jdime*nnode+inode];
			}
		}
	}
	
	return;
}