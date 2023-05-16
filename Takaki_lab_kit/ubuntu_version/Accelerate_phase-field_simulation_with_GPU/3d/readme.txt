#--------------------------------------------------------
# GPU

## Environment
- Windows11
- CPU: 12th Gen Intel(R) Core(TM) i7-12700
- GPU: NVIDIA GeForce RTX 3070 (compute capability, 8.6)
  (see, https://developer.nvidia.com/cuda-gpus )

## Get nvcc
1. sudo apt -y update
(2. sudo apt install nvidia-prime)
(3. sudo apt install cuda)
4. sudo apt -y install nvidia-cuda-toolkit
5. nvidia-smi

## Version check and adress
6. nvcc -V
7. which nvcc
  (/usr/bin/nvcc)

## Linux (ubuntu 22.04 lts) (GPU)
8. nvcc -O2 main-gpu_3d.cu write_vtk_grid_values_3D.cu -o main-gpu_3d.exe -arch=native -lm --std 'c++17'
9. ./main-gpu_3d.exe
10. (use ParaView for time_XX.vtk)

## Linux (ubuntu 22.04 lts) (CPU only)
8. nvcc -O2 main-cpu_3d.cu write_vtk_grid_values_3D.cu -o main-cpu_3d.exe -lm
9. ./main-cpu_3d.exe
10. (use ParaView for time_XX.vtk)

## cmake version (ubuntu 22.04 lts)
8. cmake -S . -B build/ -G"Unix Makefiles"
9. cmake --build build/ --target main-gpu.exe
10. cd ./build
11. ./main-gpu.exe
12. (use ParaView for time_XX.vtk)

## Reference
[T2] http://www.measlab.kit.ac.jp/nvcc.html
#--------------------------------------------------------
## Select Target Platform
- wget https://developer.download.nvidia.com/compute/cuda/repos/wsl-ubuntu/x86_64/cuda-keyring_1.0-1_all.deb
- sudo dpkg -i cuda-keyring_1.0-1_all.deb
- sudo apt-get update
- sudo apt-get -y install cuda
[S1] https://developer.nvidia.com/cuda-downloads?target_os=Linux&target_arch=x86_64&Distribution=WSL-Ubuntu&target_version=2.0&target_type=deb_network

## Environment settings
- vim ~/.bashrc
  export PATH="/usr/local/cuda/bin:$PATH"
  export LD_LIBRARY_PATH="/usr/local/cuda/lib64:$LD_LIBRARY_PATH"
[S2] https://misoji-engineer.com/archives/ubuntu20-04-cuda.html
#--------------------------------------------------------
nvcc main-cpu.cu -o main-cpu

For CPU: main-cpu.cu
For GPU (Global memory): main-gpu.cu
For GPU (Shard memory): main-shared.cu

[T1] http://www.cis.kit.ac.jp/~takaki/phase-field/ch22-23/ch22-23.html
[T2] https://www.gsic.titech.ac.jp/supercon/main/attwiki/index.php?plugin=attach&refer=SupercomputingContest2018&openfile=SuperCon2018-GPU.pdf
[T3] https://hpc-phys.kek.jp/workshop/workshop190318/aoyama_190318.pdf
[T4] http://olab.is.s.u-tokyo.ac.jp/~reiji/GPUtut2.pdf
#--------------------------------------------------------
## Memo -1-
$ nvcc -O2 main-gpu_3d.cu write_vtk_grid_values_3D.cu -o main-gpu_3d.exe -arch=native -lm --std 'c++17'
$ ./main-gpu_3d.exe
--------------------------------------------------
Device Number: 0
  Device name: NVIDIA GeForce RTX 3070
  Memory Clock Rate (KHz): 7001000
  Memory Bus Width (bits): 256
  Peak Memory Bandwidth (GB/s): 448.064000
--------------------------------------------------
istep =     0, average constration = 0.402716, total annealing time= 0.000000e+00 [s] at 623.000000 [K]
istep =  1000, average constration = 0.401688, total annealing time= 4.098417e+11 [s] at 623.000000 [K]
istep =  2000, average constration = 0.400979, total annealing time= 8.196834e+11 [s] at 623.000000 [K]
istep =  3000, average constration = 0.400993, total annealing time= 1.229525e+12 [s] at 623.000000 [K]
istep =  4000, average constration = 0.400385, total annealing time= 1.639367e+12 [s] at 623.000000 [K]
istep =  5000, average constration = 0.400023, total annealing time= 2.049208e+12 [s] at 623.000000 [K]
istep =  6000, average constration = 0.400017, total annealing time= 2.459050e+12 [s] at 623.000000 [K]
istep =  7000, average constration = 0.400019, total annealing time= 2.868892e+12 [s] at 623.000000 [K]
istep =  8000, average constration = 0.400015, total annealing time= 3.278734e+12 [s] at 623.000000 [K]
istep =  9000, average constration = 0.400000, total annealing time= 3.688575e+12 [s] at 623.000000 [K]
istep = 10000, average constration = 0.399979, total annealing time= 4.098417e+12 [s] at 623.000000 [K]
Calculation Time =    64.879 [sec]
#--------------------------------------------------------
## Memo -2-
$ nvcc -O2 main-cpu_3d.cu write_vtk_grid_values_3D.cu -o main-cpu_3d.exe -lm
$ ./main-cpu_3d.exe
--------------------------------------------------
Device Number: 0
  Device name: NVIDIA GeForce RTX 3070
  Memory Clock Rate (KHz): 7001000
  Memory Bus Width (bits): 256
  Peak Memory Bandwidth (GB/s): 448.064000
--------------------------------------------------
istep =     0, average constration = 0.402715, total annealing time= 0.000000e+00 [s] at 673.000000 [K]
istep =  1000, average constration = 0.400425, total annealing time= 6.042171e+09 [s] at 673.000000 [K]
istep =  2000, average constration = 0.401865, total annealing time= 1.208434e+10 [s] at 673.000000 [K]
istep =  3000, average constration = 0.398974, total annealing time= 1.812651e+10 [s] at 673.000000 [K]
istep =  4000, average constration = 0.400194, total annealing time= 2.416869e+10 [s] at 673.000000 [K]
istep =  5000, average constration = 0.400069, total annealing time= 3.021086e+10 [s] at 673.000000 [K]
istep =  6000, average constration = 0.400485, total annealing time= 3.625303e+10 [s] at 673.000000 [K]
istep =  7000, average constration = 0.401545, total annealing time= 4.229520e+10 [s] at 673.000000 [K]
istep =  8000, average constration = 0.401101, total annealing time= 4.833737e+10 [s] at 673.000000 [K]
istep =  9000, average constration = 0.401099, total annealing time= 5.437954e+10 [s] at 673.000000 [K]
istep = 10000, average constration = 0.400865, total annealing time= 6.042171e+10 [s] at 673.000000 [K]
Calculation Time =   696.143 [sec]
#--------------------------------------------------------