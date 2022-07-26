check ubuntu 18.04 LTS
-----------------------------------------------------------------------
1. sudo apt update
2. sudo apt -y install g++
3. sudo apt -y install gnuplot
4. g++ main-cpu.cpp -o main-cpu
5. ./main-cpu
6. gnuplot
  set pm3d map
  splot "f010.dat"
  quit
-----------------------------------------------------------------------