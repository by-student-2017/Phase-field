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
8. nvcc -O2 main-gpu.cu write_vtk_grid_values_3D.cu -o main-gpu.exe -arch=native -lm --std 'c++17'
9. ./main-gpu.exe
10. (use ParaView for time_XX.vtk)

## Linux (ubuntu 22.04 lts) (GPU, using shared memory)
8. nvcc -O2 main-shared.cu write_vtk_grid_values_3D.cu -o main-shared.exe -arch=native -lm --std 'c++17'
9. ./main-shared.exe
10. (use ParaView for time_XX.vtk)

## Linux (ubuntu 22.04 lts) (CPU only)
8. nvcc -O2 main-cpu.cu write_vtk_grid_values_3D.cu -o main-cpu.exe -lm
9. ./main-shared.exe
10. (use ParaView for time_XX.vtk)

## cmake version (ubuntu 22.04 lts) (failed)
8. cmake -S . -B build/ -G"Unix Makefiles"
9. cmake --build build/ --target main-gpu.exe
10. cd ./build
11. ./main-gpu.exe
12. (use ParaView for time_XX.vtk)

## fft version (underconstracting)
- nvcc -O2 main-gpu.cu write_vtk_grid_values_3D.cu -o main-gpu.exe -lcufft -arch=sm_86

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