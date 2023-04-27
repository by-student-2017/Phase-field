Calculation procedure

1. pfc_2d
  (get output file: final_conf.out)
2. pfc_poly_2d
  (input: final_conf.out from pfc_2d)
  (output: bi_2r_2d.inp)
  (Select a circle or vertical line defect in the code (pfc_poly_2d.c).)
3. pfc_def_2d or pfc_2d
  (input: bi_2r_2d.inp)

"3d" is the same procedure as "2d".

#----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----
M is the mobility.
F is the free energy of the system as function of phase field phi(r).
r is the position vector.
q0 is the wavenumber.
alpha*dT is the driving force.
#----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----
  x = q0 * r
  psai = sqrt(g/(lambda*q0^4)) * (phi(r) - h/(3*g))
  epsilon = 1/(lambda*q0^4) * (h/(3*g) - alpha)
where, alpha, lambda, h, q0, and g are phenomenological parameters to fit 
properties of the material of interest.
It is convenient to rewrite the Eq.7.2 in dimensionless form, 
by introducing a new set of variables as: in which all terms that are 
either constant or linearly dependent to 
the dimensionless order parameter field psai can beignored and 
the value of h is usually taken as zero.
Then, a dimensionless free energy F~ as g*lambda^-2*q0^-5 times
the original free energy, without the constant and linear parts, 
takes the form,
  (Eq.7.4)
Here, epsilon is a constant proportional to
the deviation of the temperature from the melting temperature; 
therefore, it takes a negative value.

Assumption
  constant phase: liquid phase
  triangle phase: solid phase
#----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----
