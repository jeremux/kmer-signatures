import weka.core.Instance;
import weka.core.Instances;
import weka.core.Utils;
import weka.core.converters.ArffLoader;
import weka.classifiers.Evaluation;
import weka.classifiers.bayes.NaiveBayesUpdateable;
import weka.classifiers.functions.LibSVM;
import weka.classifiers.bayes.NaiveBayes;

import java.io.File;


/**
 * This example trains NaiveBayes incrementally on data obtained
 * from the ArffLoader.
 *
 * @author FracPete (fracpete at waikato dot ac dot nz)
 */
public class IncrementalClassifier  {

	
	public static double getCorrect(String root,int id) throws Exception
	{
		double res=0;
		
		String trainPath = root + "/frequencies/learn-"+id+".arff";
		String testDataset = root + "/frequencies/toPredict-"+id+".arff";
		ArffLoader loader = new ArffLoader();
		loader.setFile(new File(trainPath));
		Instances train = loader.getStructure();
		train.setClassIndex(train.numAttributes()-1);

	
		NaiveBayes nb = new NaiveBayesUpdateable();
//		nb.setOptions(Utils.splitOptions("-K"));
		nb.buildClassifier(train);
//		int cpt = 0;
		

		Instance current;
		while ((current = loader.getNextInstance(train)) != null)
		{
//			cpt++;
//			System.out.println("cpt = "+cpt);
			nb.updateClassifier(current);
		}


		// output generated model
//		System.out.println(nb);

		// Load the test dataset
		ArffLoader loader_test = new ArffLoader();
		loader_test.setFile(new File(testDataset));

		Instances instances_test = loader_test.getDataSet();
		instances_test.setClassIndex(instances_test.numAttributes()-1);

		// Evaluate the model
		Evaluation eval = new Evaluation(instances_test);
		eval.evaluateModel(nb, instances_test);

		// Print summary
//		System.out.println(eval.toSummaryString());

		// Print confusion matrix
//		double correct = eval.correct();
//		double incorrect = eval.numInstances() - correct;
		res = (int) Math.round((eval.correct()/eval.numInstances())*100);
		
		
		return res;
	}
	public static void main(String[] args) throws Exception 
	{
		double res=0;
		int nTests = Integer.parseInt(args[1]);
		for(int i=1;i<=nTests;i++)
		{
			res+= getCorrect(args[0],i);
		}
		res /= (double) nTests;
		
		System.out.println(res);
	}
}
