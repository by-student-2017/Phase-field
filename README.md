# Phase-field

## Installation
1. (open Linux. e.g., Ubuntu or WSL(windows))
2. cd $HOME
3. git clone https://github.com/by-student-2017/Phase-field.git
4. sudo apt update
5. sudo apt -y install g++
6. sudo apt -y install paraview paraview-dev

## Tutorial 1 (basic run)
1. cd $HOME/Phase-field
2. cd Koyama_lab_nu
3. cd programs_ubuntu_version
4. cd Chapter4_Section2_three_phase_v2
5. ls
6. cat readme.txt
7. g++ three_phase_v2.cpp -o three_phase
8. ./three_phase

## copy linux to win
9. cd ../
10. ls
11. cp -r Chapter4_Section2_three_phase_v2 /mnt/c/Users/*/Desktop

## Tutorial 2 (change parameter)
12. vim parameters.txt
13. (keyboard) i
14. (change "temp= 900.0        ! Temperature [K]" to "temp= 450.0        ! Temperature [K]")
15. (keyboard) Esc
16. (keyboard) :wq
17. (keyboard) Enter
18. ./three_phase

## copy win to linux
12. (rewrite parameter of parameters.txt)
13. cp /mnt/c/Users/*/Desktop/Chapter4_Section2_three_phase_v2/parameters.txt $HOME/Phase-field/Koyama_lab_nu/programs_ubuntu_version/Chapter4_Section2_three_phase_v2/

## Note 1
- Press (keyboard) [tab] while entering the text will automatically write a continuation or search for candidates.
- "Paraview" can be used on Windows, so don't worry about it. Load the vtk file output from the calculation into "Paraview".
- Steps 7, 8, 9 and 10 only need to be done once.

## Note 2
- Please read "readme.txt" in each file for rights and citations.
- I recommend using "ubuntu 18.04 LTS" for WSL (windows10) due to image display issues.

## XLaunch (ubuntu 20.04 LTS or ubuntu 22.04 LTS on windows10)
1. sudo strip --remove-section=.note.ABI-tag /usr/lib/x86_64-linux-gnu/libQt5Core.so.5
2. echo "export DISPLAY=:0.0" >> ~/.bashrc
3. bash

## ParaView on ubuntu 20.04 LTS 
- move red line to origin (>1.0 >1.0 >1.0)

## Editor
- SAKURA: https://sakura-editor.github.io/index.en.html
- SAKURA: https://sakura-editor.github.io/ (Japanese)

## WSL (windows 11 + VcXsrv)
- see https://www5.hp-ez.com/hp/calculations/page513 (Japanese)
1. Settings -> Projection to PC -> Optional Features -> Other Windows Features
2. ☑ Windows Subsystem for Linux
3. ☑ Virtual Machine Platform
4. OK
5. (Open a command prompt with administrator privileges)
6. wsl --list --d
7. wsl --install -d Ubuntu-20.04
8. wsl --update
9. wsl --shutdown

## GPU setting (not need this setting)
- see Takaki_lab_kit/ubuntu_version/Accelerate_phase-field_simulation_with_GPU/2d/readme.txt

## Acknowledgment
- This project (modified version) is/was partially supported by the following :
  + meguREnergy Co., Ltd.
  + ATSUMITEC Co., Ltd.

## Information of CUDA vs. OpenACC
- https://proc-cpuinfo.fixstars.com/2020/06/openacc-vs-openmp-gpu-offloading/
- https://www.jsces.org/activity/journal/files/tutorial_2103.pdf
- https://www.cc.u-tokyo.ac.jp/events/lectures/209/20230630-1.pdf
- https://www.cfca.nao.ac.jp/files/20230117ymiki.pdf
