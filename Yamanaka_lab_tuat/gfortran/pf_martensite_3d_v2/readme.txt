-----------------------------------------------------------------------
ubuntu 20.04 LTS
(underconstructing)

1. sudo apt update
2. sudo apt -y install gfortran libfftw3-dev
3. sudo apt -y install paraview paraview-dev
4. gfortran pf_martensite_v2_libfftw3.f90 -I/usr/include -lfftw3 -o pf_martensite_v2
5. ./pf_martensite_v2
6. paraview
  a1. File -> Open ... -> mt_result..vtk -> OK
  a2. (click) [Apply] 
  a3. (click) Clip1
-----------------------------------------------------------------------