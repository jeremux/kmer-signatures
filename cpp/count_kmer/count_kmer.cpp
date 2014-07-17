/*
 * count_kmer.cpp
 *
 *  Created on: 10 juil. 2014
 *      Author: jeremy
 */
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <typeinfo>
#include <iomanip>
#include <getopt.h>


#include "classPattern.h"
#include "classData.h"
#include "FreqKmer.h"

#define no_argument 0
#define required_argument 1
#define optional_argument 2
#define VERSION 1.0


typedef struct {
  string listFastaPath;
  string fastaPath;
  int windowSize;
  string output;
  string kmerPath;
  Switch version;
  Switch doPrintHelp;
  Switch doTest;
} options;

options opt;
bool initFromList=false;

void init_opt()
{
  opt.listFastaPath = "null";
  opt.fastaPath = "null";
  opt.windowSize = 0;
  opt.output = "null";
  opt.kmerPath = "null";
  opt.version = No;
  opt.doPrintHelp =No;
  opt.doTest = No;
}

int getParam(int argcc, char **argvv,options *opt)
{
	init_opt();
	int c;
	while (1)
	{
		int option_index = 0;
		static struct option long_options[] = {
				{"help", no_argument, 0, 'h'},            		//-h
				{"listFasta", required_argument, 0, 'F'},       //-F
				{"fasta", required_argument, 0, 'f'},           //-f
				{"wsize", required_argument, 0, 'l'},    		//-l
				{"kmer", required_argument, 0, 'k'},            //-k
				{"output", required_argument, 0, 'o'},          //-o
				{"version", no_argument, 0, 'v'},               //-v
				{"test",no_argument, 0, 't'},        			//-t
				{0, 0, 0, 0}
		};
		c = getopt_long (argcc, argvv, "hF:f:l:k:o:vt",long_options, &option_index);
		if (c == -1)
			break;
		switch (c) {
		case 0: /* get a long option */
			if(!(strcmp(long_options[option_index].name,"help")))
			{
				opt->doPrintHelp=Yes;
			}
			if(!(strcmp(long_options[option_index].name,"test")))
			{
				opt->doTest=Yes;
			}
			if(!(strcmp(long_options[option_index].name,"version")))
			{
				opt->version=Yes;
			}
			if(!(strcmp(long_options[option_index].name,"listFasta")))
			{
				if(optarg)
				{
					opt->listFastaPath=optarg;
					initFromList = true;
				}
			}
			if(!(strcmp(long_options[option_index].name,"fasta")))
			{
				if(optarg)
					opt->fastaPath=optarg;
			}
			if(!(strcmp(long_options[option_index].name,"output")))
			{
				if(optarg)
					opt->output=optarg;
			}
			if(!(strcmp(long_options[option_index].name,"kmer")))
			{
				if(optarg)
					opt->kmerPath=optarg;
			}
			if(!(strcmp(long_options[option_index].name,"wsize")))
			{
				if(optarg)
					opt->windowSize=atoi(optarg);
			}
			break;
		case 'h':
			opt->doPrintHelp=Yes;
			break;
		case 't':
			opt->doTest=Yes;
			break;
		case 'v':
			opt->version=Yes;
			break;
		case 'F':
			if(optarg)
			{
				opt->listFastaPath=optarg;
				initFromList = true;
			}
			break;
		case 'f':
			if(optarg)
				opt->fastaPath=optarg;
			break;
		case 'o':
			if(optarg)
				opt->output=optarg;
			break;
		case 'k':
			if(optarg)
				opt->kmerPath=optarg;
			break;
		case 'l':
			if(optarg)
				opt->windowSize=atoi(optarg);
			break;
		default:
			cerr << "?? getopt returned character code " <<  c << "\n";;
		}
	}
	if (optind < argcc) {
		cerr << "non-option ARGV-elements: ";
		while (optind < argcc)
			cerr << "!" << argvv[optind++] << "!\n";

	}
	return 1;
}


void printHelp()
{
	printf("usage: count_kmer [-F listPathFasta | -f fasta_file ] -k pattern_kmer -l n -o file_out \n\n"
		 "count mers in fasta_file with specified kmer in file pattern_kmer\n"
		 "\n"
		 "  --listFasta		-F  file with a list of fasta path\n"
		"  --fasta	        -f  fasta file\n"
		 "  --wsize   		-l  length of the read (if -1 read = sequence)\n"
		 "  --kmer     		-k  file with pattern of mers to count\n"
		 "  --output   		-o  output filename\n"
		 "  --version  		-v  program version\n"
		 "  --help     		-h  print this help\n"
		 "\n"
		 " jeremy.fontaine@etudiant.univ-lille1.fr\n");

}

void printVersion()
{
	cout << "count_kerm " << VERSION << "\n";
}
void doTest3()
{
	cerr << "\n**test multi-pattern**\n";

	FreqKmer *f = new FreqKmer(4);
	string filename = "tests/test1/liste.txt";
	f->initPatterns("tests/test1/pattern.txt");
	f->initFromList(filename);

	int ligne = 7;
	int col = 20;
	double ** freq = new double *[ligne];
	for(int i=0 ; i<ligne ; i++)
	{
		freq[i] = new double[col];
		for(int j=0 ; j<col ; j++)
		{
			freq[i][j]=0.0;
		}
	}
	freq[0][0] = 1;
	freq[0][5] = 1;
	freq[1][3] = 1;
	freq[1][5] = 1;
	freq[2][3] = 1;
	freq[2][6] = 1;
	freq[3][3] = 1;
	freq[3][13] = 1;
	freq[4][12] = 1;
	freq[4][13] = 1;
	freq[5][6] = 1;
	freq[5][12] = 1;
	freq[6][6] = 1;
	freq[6][1] = 1;
	freq[0][16] = 2;
	freq[0][17] = 2;
	freq[1][16] = 1;
	freq[1][17] = 2;
	freq[1][19] = 1;
	freq[2][16] = 1;
	freq[2][17] = 1;
	freq[2][18] = 1;
	freq[2][19] = 1;
	freq[3][16] = 1;
	freq[3][17] = 1;
	freq[3][19] = 2;
	freq[4][16] = 1;
	freq[4][17] = 1;
	freq[4][19] = 2;
	freq[5][16] = 1;
	freq[5][17] = 1;
	freq[5][18] = 1;
	freq[5][19] = 1;
	freq[6][16] = 1;
	freq[6][17] = 2;
	freq[6][18] = 1;

	f->fillFreq();
//	cout << "Fin fillFreq\n";
//	cout << "Debut Impression\n";
	bool res=true;
	for(int i=0;i<ligne;i++)
	{
		for(int j=0;j<col;j++)
		{
			if(f->getFreq()[i][j]!=freq[i][j])
			{
				res=false;
				break;
			}
		}
	}
	if (res)
	{
		cerr << "Test3 OK !\n";
	}
	else
	{
		cerr << "Test3 fail...\n";
	}
	delete f;
}

void doTest2()
{
	cerr << "\n**test pattern simple**\n";
	FreqKmer *g = new FreqKmer(4);
	string filename = "tests/test1/liste.txt";
	g->initPatterns("tests/test2/pattern2.txt");
	g->initFromList(filename);

	int ligne = 7;
	int col = 16;
	double ** freq = new double *[ligne];
	for(int i=0 ; i<ligne ; i++)
	{
		freq[i] = new double[col];
		for(int j=0 ; j<col ; j++)
		{
			freq[i][j]=0.0;
		}
	}
	freq[0][0] = 1;
	freq[0][5] = 1;
	freq[1][3] = 1;
	freq[1][5] = 1;
	freq[2][3] = 1;
	freq[2][6] = 1;
	freq[3][3] = 1;
	freq[3][13] = 1;
	freq[4][12] = 1;
	freq[4][13] = 1;
	freq[5][6] = 1;
	freq[5][12] = 1;
	freq[6][6] = 1;
	freq[6][1] = 1;

	g->fillFreq();
//	cout << "Fin fillFreq\n";
//	cout << "Debut Impression\n";
	bool res=true;
	for(int i=0;i<ligne;i++)
	{
		for(int j=0;j<col;j++)
		{
			if(g->getFreq()[i][j]!=freq[i][j])
			{
				res=false;
				break;
			}
		}
	}
	if (res)
	{
		cerr << "Test2 OK !\n";
	}
	else
	{
		cerr << "Test2 fail...\n";
	}
	delete g;
}

void doTest1()
{
	cerr << "\n**test dimension tableau**\n";
	FreqKmer *f = new FreqKmer(6);
	string filename = "tests/test3/test.fasta";

	f->initPatterns("tests/test3/3patterns.txt");
	f->initFromFasta(filename);
	if(f->getNCol()!=84 || f->getNLigne()!=107)
	{
		cerr << "Test1 fail...\n";
	}
	else
	{
		cerr << "Test1 OK !\n";
	}

	delete f;
}

void callTest(int i)
{
	switch (i) {
		case 1:
			doTest1();
			break;
		case 2:
			doTest2();
			break;
		case 3:
			doTest3();
			break;
		default:
			break;
	}
}



void executeTests(int k)
{
	cerr << "===Début des tests===\n";
	for(int i=1; i <= k ; i++)
	{
		callTest(i);
	}
	cerr << "\n===Fin des tests===\n";
}


int main(int argc, char **argv) {

	getParam(argc,argv,&opt);
	bool error = false;
	if (opt.doPrintHelp==Yes)
	{
	  printHelp();
	  exit(0);
	}
	if (opt.doTest)
	{
		executeTests(3);
		exit(0);
	}

	if (opt.listFastaPath=="null" && opt.listFastaPath=="null")
	{
		cerr << "Need a fasta file (-f file.fasta) or a list of fasta path (-F listFasta_file) \n";
		error = true;
	}

	if (opt.kmerPath=="null")
	{
		cerr << "Need a file with kmer pattern (-k kmer_file)\n";
		error = true;
	}

	if (opt.output=="null")
	{
		cerr << "Need a filename for the output (-o outputFile)\n";
		error = true;
	}
	if (opt.windowSize == 0)
	{
		cerr << "Need the size of the read (-l read_size)\n";
		error = true;
	}

	if (error)
		exit(1);
	FreqKmer *f = new FreqKmer(opt.windowSize);

	if(initFromList)
	{
		f->initFromList(opt.listFastaPath);
	}
	else
	{
		f->initFromFasta(opt.fastaPath);
	}
	f->initPatterns(opt.kmerPath);

	f->fillFreq();
	f->imprimeCSV(opt.output);

	return 0;
}





