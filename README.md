# Phase-field

## Installation
1. (open Linux)
2. cd $HOME
3. git clone https://github.com/by-student-2017/Phase-field.git

## Tutorial 1 (basic run)
1. cd $HOME/Phase-field
2. cd Koyama_lab_nu
3. cd programs_ubuntu_version
4. cd Chapter4_Section2_three_phase_v2
5. ls
6. cat readme.txt
7. sudo apt update
8. sudo apt -y install g++
9. sudo apt -y install paraview paraview-dev
10. g++ three_phase_v2.cpp -o three_phase
11. ./three_phase

## copy linux to win
12. cd ../
13. ls
14. cp -r Chapter4_Section2_three_phase_v2 /mnt/c/Users/*/Desktop

## Tutorial 2 (change parameter)
15. vim parameters.txt
16. (keyboard) i
17. (change "temp= 900.0        ! Temperature [K]" to "temp= 450.0        ! Temperature [K]")
18. (keyboard) Esc
19. (keyboard) :wq
20. (keyboard) Enter
21. ./three_phase

## copy win to linux
15. (rewrite parameter of parameters.txt)
16. cp /mnt/c/Users/*/Desktop/Chapter4_Section2_three_phase_v2/parameters.txt $HOME/Phase-field/Koyama_lab_nu/programs_ubuntu_version/Chapter4_Section2_three_phase_v2/

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

##Acknowledgment ######################################
- This project (modified version) is/was partially supported by the following :
  + meguREnergy Co., Ltd.
  + ATSUMITEC Co., Ltd.

## Information of CUDA vs. OpenACC
- https://proc-cpuinfo.fixstars.com/2020/06/openacc-vs-openmp-gpu-offloading/
- https://www.jsces.org/activity/journal/files/tutorial_2103.pdf
- https://www.cc.u-tokyo.ac.jp/events/lectures/209/20230630-1.pdf
- https://www.cfca.nao.ac.jp/files/20230117ymiki.pdf
