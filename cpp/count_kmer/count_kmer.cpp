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
#include "test.h"

#define no_argument 0
#define required_argument 1
#define optional_argument 2
#define VERSION 1.0

#define NB_TEST 10

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
        executeTests(NB_TEST);
//		FreqKmer *f = new FreqKmer(5,2);
//		string filename = "test.fasta";
//		f->initPatterns("pattern.txt");
//		f->initFromFasta(filename);
//
//		f->fillFreq();
//		f->imprimeCSV("tata.csv");
//		doTest10();

		exit(0);
	}

	if (opt.listFastaPath=="null" && opt.fastaPath=="null")
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
		f->initDataFromListFastaPath(opt.listFastaPath);
	}
	else
	{
		cerr << "Debut init fasta \n";
		f->initFromFasta(opt.fastaPath);
		cerr << "Fin init fasta\n";
	}
	cerr << "Debut init patterns\n";
	f->initPatterns(opt.kmerPath);
	cerr << "Fin init patterns\n";

	cerr << "Debut fill \n";
	f->fillFreq();
	cerr << "Fin fill \n";
	f->imprimeCSV(opt.output);

	return 0;
}





