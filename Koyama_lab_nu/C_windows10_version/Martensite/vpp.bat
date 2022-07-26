set path_1=%path%
set VS6=..\Visual_Studio
set LIB=%VS6%\Vc98\lib;
path %VS6%\Vc98\Bin;%VS6%\Vc98\Include;%VS6%\Common\MsDev98\Bin
CL /I%VS6%\slswin /I%VS6%\Vc98\include %1 %2 %3 %4 %5 %6 %7 %8 %9 /link kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib
set LIB=
path ;
path %path_1%

