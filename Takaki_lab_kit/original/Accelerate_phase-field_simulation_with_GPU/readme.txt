
nvcc main-cpu.cu -o main-cpu

For CPU: main-cpu.cu
For GPU (Global memory): main-gpu.cu
For GPU (Shard memory): main-shared.cu

[T1] http://www.cis.kit.ac.jp/~takaki/phase-field/ch22-23/ch22-23.html
#--------------------------------------------------------
# GPU

## Environment
- Windows11
- CPU: 12th Gen Intel(R) Core(TM) i7-12700
- GPU: NVIDIA GeForce RTX 3070 (compute capability, 8.6)
  (see, https://developer.nvidia.com/cuda-gpus )

## Get nvcc
- sudo apt -y update
- sudo apt -y install nvidia-cuda-toolkit
- nvidia-smi

## Version check and adress
- nvcc -V
- which nvcc
  (/usr/bin/nvcc)

## Linux (ubuntu 22.04 lts)
- nvcc -O2 main-gpu.cu -lcutil -o main-gpu.exe
- ./main-gpu.exe

## fft version
- nvcc -O2 main-gpu.cu -o main-gpu.exe -lcufft

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