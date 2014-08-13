javac -classpath $CLASSPATH:../../weka-3-6-11/weka.jar:../../weka-3-6-11/libsvm.jar IncrementalClassifier.java
java -classpath $CLASSPATH:../../weka-3-6-11/weka.jar:../../weka-3-6-11/libsvm.jar IncrementalClassifier $1 $2
