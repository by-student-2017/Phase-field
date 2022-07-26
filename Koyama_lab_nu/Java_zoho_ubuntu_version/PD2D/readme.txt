check ubuntu 18.04 LTS
-----------------------------------------------------------------------
1. sudo apt update
2. sudo apt -y install default-jre
3. java PD2D_001
-----------------------------------------------------------------------
- other method -1-
1. sudo apt update
2. sudo apt -y install default-jre
3. sudo apt -y install default-jdk
4. sudo apt -y install nkf
5. nkf -g PD2D_001.java
6. nkf -w --overwite FePt2D_001.java
7. javac PD2D_001.java -Xlint
8. rm -f *.class
9. java PD2D_001

1. nkf -w --overwite plot_c.java
2. javac plot_c.java -Xlint

1. nkf -w --overwite bplot_c.java
2. javac bplot_c.java -Xlint
-----------------------------------------------------------------------
- other method -2-
1. sudo apt update
2. sudo apt -y install default-jre
3. sudo apt -y install default-jdk
4. iconv -f SJIS -t UTF8 PD2D_001.java > PD2D_001_new.java
5. rm -f PD2D_001.java *.class 
6. mv PD2D_001_new.java PD2D_001.java
7. javac PD2D_001.java -Xlint
8. java PD2D_001

1. iconv -f SJIS -t UTF8 plot_c.java > plot_c_new.java
2. rm -f plot_c.java
3. mv plot_c_new.java plot_c.java
4. javac plot_c.java -Xlint

1. iconv -f SJIS -t UTF8 bplot_c.java > bplot_c_new.java
2. rm -f bplot_c.java
3. mv bplot_c_new.java bplot_c.java
4. javac bplot_c.java -Xlint
-----------------------------------------------------------------------