check ubuntu 18.04 LTS
-----------------------------------------------------------------------
1. sudo apt update
2. sudo apt -y install default-jre
3. java java F_curve_FePt
-----------------------------------------------------------------------
- other method -1-
1. sudo apt update
2. sudo apt -y install default-jre
3. sudo apt -y install default-jdk
4. sudo apt -y install nkf
5. nkf -g F_curve_FePt.java
6. nkf -w --overwite FePt2D_001.java
7. javac F_curve_FePt.java -Xlint
8. rm -f *.class
9. java F_curve_FePt
-----------------------------------------------------------------------
- other method -2-
1. sudo apt update
2. sudo apt -y install default-jre
3. sudo apt -y install default-jdk
4. iconv -f SJIS -t UTF8 F_curve_FePt.java > F_curve_FePt_new.java
5. rm -f F_curve_FePt.java *.class 
6. mv F_curve_FePt_new.java F_curve_FePt.java
7. javac F_curve_FePt.java -Xlint
8. java F_curve_FePt
-----------------------------------------------------------------------