-----------------------------------------------------------------------
ubuntu 18.04 LTS
(development version)
(Note: Due to lack of time, it could not be fully developed. 
  If there are any mistakes, please point them out and correct them.)

1. sudo apt update
2. sudo apt -y install g++
3. sudo apt -y install paraview paraview-dev
4. g++ dend_4fold_v2.cpp -o dend_4fold
5. ./dend_4fold
6. paraview
  a1. File -> Open ... -> sp_result..vtk -> OK
  a2. (click) [Apply] 
  a3. (click "play") |>

[1] T. Takaki, and A. Yamanaka., "Phase field method", Yokendo (Japanese)
  ISBN-10 : 4842504927
  ISBN-13 : 978-4842504926
  See Chapter 4 for more information on expressions.
-----------------------------------------------------------------------
