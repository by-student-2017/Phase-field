ubuntu 18.04 LTS
-----------------------------------------------------------------------
1. sudo apt update
2. sudo apt -y install gfortran gnuplot
3. gfortran ch14_mpf_gnuplot.f -o mpf
4. ./mpf
5. gnuplot
  set pm3d map
  splot "out000500.dat" u 1:2:3
  quit
-----------------------------------------------------------------------