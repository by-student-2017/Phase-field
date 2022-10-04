
programs_Java.zip, programs_C.zip and 3D_programs.zip
Copyright (C) 2011, Nagoya University, All Rights Reserved.

programs_Java_zoho.zip and PF_programs_Python.zip
Copyright (C) 2019, Nagoya University, All Rights Reserved.

Reference
[1] https://www.material.nagoya-u.ac.jp/PFM/Phase-Field_Modeling.htm

google_colab is google colaboratory version from PF_programs_Python.zip

----- dimensionless -----
Energy: 1/RT [1/(J/mol)]
Length: 1/b1 = 1/(L/N) = 1/([m]/(number of mesh))
Time: 1/(b1*b1/D), diffusion coefficient D [m2/s]
Time: Ms*RT, Ms [mol/s/J]=[(1/s)*1/(J/mol)], (mobility Mc=Ms=1)
Mc*RT/D [dimensionless]
Concentration gradient energy, kappa: 1/(b1*b1*RT)

Einstein relations: D = M0*RT, M0 [(m2/s)/(J/mol)]
-----------------------------------------------------------------------
- In the difference lattice, the wavenumber is the coordinate number of the array as it is, so if the difference division is 2^7 = 128, k = 0,1,2,... It becomes 127.
- In periodic boundary conditions, if the difference division is 2^7 = 128, it is the same up to k = 0, 1, 2,...,63, but from 64 it must be (k-128) (that is, -64, -63, ..., -1).

- If the function f(r) after the inverse Fourier transform (IFFT) is a real function, then the real part of f(k) is an even function, and the imaginary part of f(k) is an odd function. This means f(-k)=f*(k). * means complex conjugate).
- In the calculation of the displacement field uj(k), the imaginary number i is multiplied, so that -Zj(k)*phi(k) is assigned to the real part and Zj(k)*phi(k) is assigned to the imaginary part.
- In the actual program, in order to prioritize the readability of the program based on the symmetry of the mathematical formula, -Zj(k)*[real(phi(k))+imag(phi(k))]=-Zj(k)*phi(k) is assigned to the real part and Zj(k)*[real(phi(k))+imag(phi(k))]=Zj(k)*phi(k) is assigned to the imaginary part. The extra terms -Zj(k)*real(phi(k)) of the real part (which is an odd function) and the imaginary Zj(k)*imag(phi(k)) (which is an even function) disappear during the inverse Fourier transform (IFFT).
-----------------------------------------------------------------------