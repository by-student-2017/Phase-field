# Phase-field

## Features
- This github provides examples for not only "libfftw3" but also "cuFFT". By writing it in C, we were able to expand it to use "cuFFT". It is characterized by not only two -dimensional but also a 3D version of code.
- The C code is constructed according to the formulas in Dr. Biner's textbook (https://link.springer.com/book/10.1007/978-3-319-41196-5). There are also several settings for Paraview that are not described in the textbook, so please refer to them. The Paraview settings for the progression of destruction in "Phase-Field Crystal (PFC) method" are particularly worthy of mention.
- The Phase-field method code shown here is written exactly as in the textbook, so it is easier to understand than the DVXalpha method SCAT, FLAPW method as represented by WIEN2k, and KKR method codes such as Akai-KKR. The Phase-field method code is also suitable for learning coding techniques

- Biner_lab_inl is a C-language, GPU-compatible version of the textbook code by Dr. Biner.
- Koyama_lab_nu is a textbook code by Dr. Koyama that reads parameters from the text. It is also characterized by the output for "Paraview".

- dynamic memory is dynamic memory allocation that allows calculations of 4 GB or more (i.e. allows calculations of large systems)
- static memory is static memory allocation that allows calculations of less than 4 GB (only small systems can be calculated, but the calculation speed is slightly faster. In addition, the code is a little simpler, so it is suitable for initial learning. If you want to flexibly change the scale of calculations without frequently repeating compilation, we recommend using "dynamic memory".)

- For "python" code examples, see Dr. Yamanaka's https://github.com/Yamanaka-Lab-TUAT or Dr. Koyama's https://www.material.nagoya-u.ac.jp/PFM/Phase-Field_Modeling.htm (e.g., https://www.material.nagoya-u.ac.jp/PFM/zairyougaku_zoho_koyama.htm) (Japanese)

- I plan to continue using C language in the future. Python has many advantages, but it often causes confusion due to specification changes. Readers who are considering long -term use need to consider the use of different languages, such as C and Fortran and Octave.

- Undeveloped: 3D dendrite growth, effects of oxygen defects in PbTiO3, crystal growth dependence on plane orientation, etc.

## Installation
1. (open Linux. e.g., Ubuntu or WSL(windows))
2. sudo apt update
3. sudo apt -y install g++
4. sudo apt -y install gcc build-essential libfftw3-dev
5. cd $HOME
6. git clone https://github.com/by-student-2017/Phase-field.git
7. ls

![Installation](https://github.com/by-student-2017/Phase-field/blob/main/Fig/installation.png)
Fig. Building a computing environment using WSL2 (ubuntu 22.04 LTS on Windows). The content displayed will vary depending on the user's environment, but as long as no errors are displayed, that's fine.

## Tutorial 1 (basic run)
1. cd $HOME/Phase-field
2. cd Koyama_lab_nu
3. cd programs_ubuntu_version
4. cd Chapter4_Section2_three_phase_v2
5. ls
6. cat readme.txt
7. g++ three_phase_v2.cpp -o three_phase
8. ./three_phase

![Tutorial1](https://github.com/by-student-2017/Phase-field/blob/main/Fig/tutorial_1_basic_run.png)
Fig. Go to the directory you want to calculate, compile it in C++, and run it.　If you use the "Tab" key while entering a command, it will write the continuation or candidates. If you want to reuse previous input, the up and down keys are also useful, so it's good to remember them. The above is done using such a function.

## Copy linux to windows
9. cd ../
10. ls
11. cp -r Chapter4_Section2_three_phase_v2 /mnt/c/Users/*/Desktop
- "*" is your PC name. You can search with the [Tab] key (on your keyboard).

![Copy_linux_to_windows](https://github.com/by-student-2017/Phase-field/blob/main/Fig/copy_linux_to_windows.png)
Fig. The processing after the calculation is performed. In this example, it finished after 49 iterations. As long as the necessary data was calculated before Aborted (core dumped), there is no need to worry about it for now. Press the [Tab] key in "cp -r Chapter4_Section2_three_phase_v2 /mnt/c/Users/" to display candidates. In my environment, I have a desktop in "manab", so I save the data to the desktop as shown in the figure.

## Paraview
12. (Open "Paraview" on windows)
13. File -> Open ... -> sp_result..vtk -> OK
14. (click) [Apply]
15. (click "play") |>

![Paraview_1](https://github.com/by-student-2017/Phase-field/blob/main/Fig/paraview_1.png)
Fig. A sequence of operations for "Paraview".

![Paraview_2](https://github.com/by-student-2017/Phase-field/blob/main/Fig/paraview_2.png)
Fig. A sequence of operations for "Paraview". The data with the extension "vtk" is the data to be loaded by "Paraview". You can also drag and drop it, but this way the series of numbered data will be loaded properly.

![Paraview_3](https://github.com/by-student-2017/Phase-field/blob/main/Fig/paraview_3.png)
Fig. A sequence of operations for "Paraview". In the last image, click the "<-->" in the second or third row at the top left (near the left side of "concentration_A") to automatically change the color gradation range to an appropriate value. The figure shows "Time: 7" (=sp_result000007.vtk).

## Tutorial 2_1 (change parameter using "Copy windows to linux" method)
12. (rewrite parameter of parameters.txt) (e.g., change "temp= 900.0        ! Temperature [K]" to "temp= 450.0        ! Temperature [K]")
13. cp /mnt/c/Users/*/Desktop/Chapter4_Section2_three_phase_v2/parameters.txt $HOME/Phase-field/Koyama_lab_nu/programs_ubuntu_version/Chapter4_Section2_three_phase_v2/

## Tutorial 2_2 (change parameter on vim)
12. vim parameters.txt
13. (keyboard) i
14. (change "temp= 900.0        ! Temperature [K]" to "temp= 450.0        ! Temperature [K]")
15. (keyboard) Esc
16. (keyboard) :wq
17. (keyboard) Enter
18. ./three_phase
- Note: It is easier to rewrite using "gedit", so it is a good idea to make use of "VcXsrv" or similar.

## Q&A (What should I do with the parameters?)
- The parameters can be obtained from the TDB file used in the phase diagram calculation (CALPHAD). It is difficult to understand the contents of the TDB file, but the results are worth the effort. https://www5.hp-ez.com/hp/calculations/page516 provides a link to the paper that describes the parameters, so please use it as a reference when reading the TDB file. There are many papers in Japanese. This was not intentional, and the main reason is that I was unable to find many English papers in my environment, even though I looked for English papers.
- For information on TDB files, see "https://github.com/by-student-2017/pycalphad-v.0.10.3-examples" and "https://github.com/by-student-2017/OpenCALPHAD-v6-examples".
- If you need to use information on phases that are not included in the TDB file, you will need to obtain parameters using first-principles calculations, molecular dynamics, etc. These should be used with the understanding that they aid in qualitative or semi-quantitative understanding rather than faithfully reproducing experimental results.
- Recently, a method has been developed in which machine learning (particle methods such as Bayesian optimization) is used to find appropriate parameters while comparing images obtained from experiments with images from the phase-field method. Since it is difficult to achieve a high degree of fit with the phase-field method, like the Rietveld method with XRD, for various reasons, it is necessary to approach the method with the understanding that it is sufficient to reproduce the experiment to a certain extent.
- The values ​​of representative parameters and methods for calculating approximate values ​​are described in "ミクロ組織の熱力学" (Thermodynamics of Microstructures) (Japan Institute of Metals)(ISBN-10: 488903028X, ISBN-13: 978-4889030280)(Japanese). The subjects and their values ​​are shown in diagrams, so if your local library has it, it's a good idea to take a look at it.

## Note 1
- Press (keyboard) [tab] while entering the text will automatically write a continuation or search for candidates.
- "Paraview" can be used on Windows, so don't worry about it. Load the vtk file output from the calculation into "Paraview".

## Note 2
- Please read "readme.txt" in each file for rights and citations.
- I recommend using "ubuntu 18.04 LTS" for WSL (windows10) due to image display issues.

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
- Note: If you have "VcXsrv" or similar available, you can edit code and parameters with "gedit" instead of "vim", which is very convenient. If you find it troublesome to use "copy -r" just to edit text, it is a good idea to have "VcXsrv" or similar available (e.g., XLaunch, etc).

## XLaunch (ubuntu 20.04 LTS or ubuntu 22.04 LTS on windows10)
1. sudo strip --remove-section=.note.ABI-tag /usr/lib/x86_64-linux-gnu/libQt5Core.so.5
2. echo "export DISPLAY=:0.0" >> ~/.bashrc
3. bash

## OS
- This github code can be used in various environments such as Windows and Linux. Windows and Linux each have their own advantages and disadvantages, so please use it in the environment that is easiest for you. It is also worth considering using WSL or cygwin.

- ## Problems with line break symbols depending on the OS
- The line break symbol may differ depending on the OS, so it may not work. To change the line break symbol, do the following:
- Convert CRLF (Windows) => LF (Linux): sed -i 's/\r//g' *.txt
- Convert LF (Linux) => CRLF (Windows): sed -i 's/$/\r/g' *.txt

## Cygwin
- see https://www5.hp-ez.com/hp/calculations/page424 

## Paraview
- https://www.paraview.org/download/

## ParaView on ubuntu 20.04 LTS 
- move red line to origin (>1.0 >1.0 >1.0)

## Editor
- SAKURA: https://sakura-editor.github.io/index.en.html (free)
- SAKURA: https://sakura-editor.github.io/ (Japanese) (free)

## Phase-Field Crystal (PFC) method
- https://github.com/by-student-2017/Phase-field/tree/main/Biner_lab_inl/chapter7/dynamic_memory
- "PFC" performs calculations in the order shown below. The previous calculation serves as the input file for the next calculation.
1. pfc_2d (get output file: final_conf.out)
2. pfc_poly_2d (input: final_conf.out from pfc_2d) (output: bi_2r_2d.inp) (Select a circle or vertical line defect in the code (pfc_poly_2d.c).)
3. pfc_def_2d or pfc_2d (input: bi_2r_2d.inp)
- "3d" is the same procedure as "2d".
- Note: The "PFC" code is available on this github as well as on "https://github.com/eimrek/phase-field-crystal-mpi". However, this is quite difficult code. Before trying to improve this difficult code, I recommend that you deepen your understanding by using the code on github here and textbooks such as those by Dr. Biner.

## GPU setting (not need this setting)
- see Takaki_lab_kit/ubuntu_version/Accelerate_phase-field_simulation_with_GPU/2d/readme.txt
- Takaki_lab_kit/ubuntu_version/Accelerate_phase-field_simulation_with_GPU
- "cufft" in Biner_lab_inl
- Note: Although it depends on the calculation conditions and the system, a GPU is about 10 times faster than parallelization using a general high-end CPU. If you have the time, it is a good idea to use it. In normal use, this difference is often wiped out when you are taking a nap, so it is a good idea to consider using a GPU when you are experiencing problems such as a calculation not completing on the CPU.

## PC specs used for test
+ OS: Microsoft Windows 11 Home 64 bit
+ BIOS: 1.14.0
+ CPU： 12th Gen Intel(R) Core(TM) i7-12700
+ Base Board：0R6PCT (A01)
+ Memory：32 GB
+ GPU: NVIDIA GeForce RTX3070
+ WSL2: VERSION="22.04.1 LTS (Jammy Jellyfish)"
+ Paraview 5.11.2
+ Python 3.10.12

## Appendix: from Windows to Linux
1. Code => Download ZIP
2. (unpack) Phase-field-main.zip
3. Change the file name in the folder to "Phase-field" and place it on the desktop.
4. (open WSL, e.g., ubuntu)
5. cd $HOME
6. cp -r  /mnt/c/Users/*/Desktop/Phase-field ./
7. ls
- "/mnt/c" = C drive
- The asterisk (*) is "something appropriate (here PC or username)
- "-r" is an option that means to copy the contents of the directory (=folder)

## Acknowledgment
- This project (modified version) is/was partially supported by the following :
  + meguREnergy Co., Ltd.
  + ATSUMITEC Co., Ltd.

## Information of CUDA vs. OpenACC
- https://proc-cpuinfo.fixstars.com/2020/06/openacc-vs-openmp-gpu-offloading/
- https://www.jsces.org/activity/journal/files/tutorial_2103.pdf
- https://www.cc.u-tokyo.ac.jp/events/lectures/209/20230630-1.pdf
- https://www.cfca.nao.ac.jp/files/20230117ymiki.pdf
