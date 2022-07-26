ubuntu 18.04 LTS
(windows needs XLaunch)
Note: Comments in the source code are in Japanese.
-----------------------------------------------------------------------
1. sudo apt update
2. sudo apt -y install default-jre
----------------------------------
3. java FePt2D_001
----------------------------------
3. java FePt2D_A001
----------------------------------
-----------------------------------------------------------------------
- other method -1-
1. sudo apt update
2. sudo apt -y install default-jre
3. sudo apt -y install default-jdk
4. sudo apt -y install nkf
5. nkf -g FePt2D_001.java
6. nkf -w --overwite FePt2D_001.java
7. javac FePt2D_001.java -Xlint
8. rm -f *.class
9. java FePt2D_001

1. nkf -w --overwite plot_s12.java
2. javac plot_s12.java -Xlint

1. nkf -w --overwite bplot_s12.java
2. javac bplot_s12.java -Xlint
-----------------------------------------------------------------------
- other method -2-
1. sudo apt update
2. sudo apt -y install default-jre
3. sudo apt -y install default-jdk
4. iconv -f SJIS -t UTF8 FePt2D_001.java > FePt2D_001_new.java
5. rm -f FePt2D_001.java *.class 
6. mv FePt2D_001_new.java FePt2D_001.java
7. javac FePt2D_001.java -Xlint
8. java FePt2D_001

1. iconv -f SJIS -t UTF8 plot_s12.java > plot_s12_new.java
2. rm -f plot_s12.java
3. mv plot_s12_new.java plot_s12.java
4. javac plot_s12.java -Xlint

1. iconv -f SJIS -t UTF8 bplot_s12.java > bplot_s12_new.java
2. rm -f bplot_s12.java
3. mv bplot_s12_new.java bplot_s12.java
4. javac bplot_s12.java -Xlint
-----------------------------------------------------------------------