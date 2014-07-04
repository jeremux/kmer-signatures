import weka.core.Instance;
import weka.core.Instances;
import weka.core.converters.ArffLoader;
import weka.classifiers.Evaluation;
import weka.classifiers.bayes.NaiveBayesUpdateable;
import weka.classifiers.meta.FilteredClassifier;
import weka.classifiers.bayes.NaiveBayes;
import weka.filters.unsupervised.attribute.StringToWordVector;

import java.io.File;
import java.io.FileReader;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.Random;

/**
 * This example trains NaiveBayes incrementally on data obtained
 * from the ArffLoader.
 *
 * @author FracPete (fracpete at waikato dot ac dot nz)
 */
public class IncrementalClassifier {

  /**
   * Expects an ARFF file as first argument (class attribute is assumed
   * to be the last attribute).
   *
   * @param args        the commandline arguments
   * @throws Exception  if something goes wrong
   */
  public static void main(String[] args) throws Exception {
    // load data
//    ArffLoader loader = new ArffLoader();
//    loader.setFile(new File("/home/jeremy/taxonomic-classification/create_db/Alveolata__33630/learn/learning_S50.arff"));
	  ArffLoader loader = new ArffLoader();
	    loader.setFile(new File("/home/jeremy/taxonomic-classification/create_db/Alveolata__33630/learn/learning_S50.arff"));
	  System.out.println("TOTO");
	  Instances train = loader.getStructure();
	  System.out.println("TATA");
	  train.setClassIndex(train.numAttributes()-1);
	  
	  String testDataset = "/home/jeremy/taxonomic-classification/create_db/Alveolata__33630/learn/learning_S50.arff";

	  NaiveBayes nb = new NaiveBayesUpdateable();
	  nb.buildClassifier(train);
	  int cpt = 0;

	    Instance current;
	    while ((current = loader.getNextInstance(train)) != null)
	    {
	    	cpt++;
	    	System.out.println("cpt = "+cpt);
	    	nb.updateClassifier(current);
	    }
	      

	    // output generated model
	    System.out.println(nb);
	    
	    // Print statistics
	    System.out.println(nb);

	    // Load the test dataset
	    ArffLoader loader_test = new ArffLoader();
	    loader_test.setFile(new File(testDataset));

	    Instances instances_test = loader_test.getDataSet();
	    instances_test.setClassIndex(16);

	    // Evaluate the model
	    Evaluation eval = new Evaluation(instances_test);
	    eval.evaluateModel(nb, instances_test);

	    // Print summary
	    System.out.println(eval.toSummaryString());

	    // Print confusion matrix
	    System.out.println("Confusion Matrix:");

	    double[][] confMatrix = eval.confusionMatrix();

	    for(int i=0; i<confMatrix.length;i++) {
	    for (int j=0;j<confMatrix[i].length;j++) {
	    System.out.print(confMatrix[i][j] + " ");
	    }
	    System.out.println();
	    }

	    }
  }
