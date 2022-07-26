check ubuntu 18.04 LTS
-----------------------------------------------------------------------
1. sudo apt update
2. sudo apt -y install default-jre
3. java PD1D_001
-----------------------------------------------------------------------
- other method -1-
1. sudo apt update
2. sudo apt -y install default-jre
3. sudo apt -y install default-jdk
4. sudo apt -y install nkf
5. nkf -g PD1D_001.java
6. nkf -w --overwite FePt2D_001.java
7. javac PD1D_001.java -Xlint
8. rm -f *.class
9. java PD1D_001

1. nkf -w --overwite plot1D_c.java
2. javac plot1D_c.java -Xlint

1. nkf -w --overwite select_data1D_c.java
2. javac select_data1D_c.java -Xlint
-----------------------------------------------------------------------
- other method -2-
1. sudo apt update
2. sudo apt -y install default-jre
3. sudo apt -y install default-jdk
4. iconv -f SJIS -t UTF8 PD1D_001.java > PD1D_001_new.java
5. rm -f PD1D_001.java *.class 
6. mv PD1D_001_new.java PD1D_001.java
7. javac PD1D_001.java -Xlint
8. java PD1D_001

1. iconv -f SJIS -t UTF8 plot1D_c.java > plot1D_c_new.java
2. rm -f plot1D_c.java
3. mv plot1D_c_new.java plot1D_c.java
4. javac plot1D_c.java -Xlint

1. iconv -f SJIS -t UTF8 select_data1D_c.java > select_data1D_c_new.java
2. rm -f select_data1D_c.java
3. mv select_data1D_c_new.java select_data1D_c.java
4. javac select_data1D_c.java -Xlint
-----------------------------------------------------------------------