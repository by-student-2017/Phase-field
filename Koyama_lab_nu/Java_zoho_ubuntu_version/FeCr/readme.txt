check ubuntu 18.04 LTS
-----------------------------------------------------------------------
1. sudo apt update
2. sudo apt -y install default-jre
3. java FeCr_PD_2D_001
-----------------------------------------------------------------------
- other method -1-
1. sudo apt update
2. sudo apt -y install default-jre
3. sudo apt -y install default-jdk
4. sudo apt -y install nkf
5. nkf -g FeCr_PD_2D_001.java
6. nkf -w --overwite FePt2D_001.java
7. javac FeCr_PD_2D_001.java -Xlint
8. rm -f *.class
9. java FeCr_PD_2D_001

1. nkf -w --overwite plot_s12.java
2. javac plot_s12.java -Xlint

1. nkf -w --overwite bplot_s12.java
2. javac bplot_s12.java -Xlint
-----------------------------------------------------------------------
- other method -2-
1. sudo apt update
2. sudo apt -y install default-jre
3. sudo apt -y install default-jdk
4. iconv -f SJIS -t UTF8 FeCr_PD_2D_001.java > FeCr_PD_2D_001_new.java
5. rm -f FeCr_PD_2D_001.java *.class 
6. mv FeCr_PD_2D_001_new.java FeCr_PD_2D_001.java
7. javac FeCr_PD_2D_001.java -Xlint
8. java FeCr_PD_2D_001

1. iconv -f SJIS -t UTF8 plot_s12.java > plot_s12_new.java
2. rm -f plot_s12.java
3. mv plot_s12_new.java plot_s12.java
4. javac plot_s12.java -Xlint

1. iconv -f SJIS -t UTF8 bplot_s12.java > bplot_s12_new.java
2. rm -f bplot_s12.java
3. mv bplot_s12_new.java bplot_s12.java
4. javac bplot_s12.java -Xlint
-----------------------------------------------------------------------