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
#include "testPart2.h"
#include "test_intramacro.h"

#define no_argument 0
#define required_argument 1
#define optional_argument 2
#define VERSION 1.0

#define NB_TEST 19

typedef struct {
    string listFastaPath;
    string fastaPath;
    int windowSize;
    string output;
    string kmerPath;
    Switch version;
    Switch doPrintHelp;
    Switch doTest;
    Switch doTestIntra;
    bool noData;
    string key;
    string root;
    int jump;
} options;

options opt;
bool initFromList=false;
bool racineBool = false;

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
    opt.doTestIntra = No;
    opt.noData = false;
    opt.key = "genomes";
    opt.root="null";
    opt.jump=0;
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
	    {"intra",no_argument, 0, 'T'},        			//-T
	    {"noData",no_argument,0,'d'},					//-d
	    {"key",no_argument,0,'K'},					    //-K
	    {"root",no_argument,0,'r'},					    //-r
	    {"jump",no_argument,0,'j'},					    //-j
	    {0, 0, 0, 0}
	};
	c = getopt_long (argcc, argvv, "hF:f:l:k:o:vtdTKr:j:",long_options, &option_index);
	if (c == -1)
	    break;
	switch (c) {
	case 0: /* get a long option */
	    if(!(strcmp(long_options[option_index].name,"help")))
	    {
		opt->doPrintHelp=Yes;
	    }
	    if(!(strcmp(long_options[option_index].name,"noData")))
	    {
		opt->noData = true;
	    }
	    if(!(strcmp(long_options[option_index].name,"test")))
	    {
		opt->doTest=Yes;
	    }
	    if(!(strcmp(long_options[option_index].name,"intra")))
	    {
	    	opt->doTestIntra=Yes;
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
	    if(!(strcmp(long_options[option_index].name,"root")))
		{
		if(optarg)
		{
			opt->root=optarg;
			racineBool=true;
			initFromList=false;
		}
		}
	    if(!(strcmp(long_options[option_index].name,"key")))
		{
		if(optarg)
		{
			opt->key=optarg;
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
	    if(!(strcmp(long_options[option_index].name,"jump")))
		{
		if(optarg)
			opt->jump=atoi(optarg);
		}
	    break;
	case 'h':
	    opt->doPrintHelp=Yes;
	    break;
	case 't':
	    opt->doTest=Yes;
	    break;
	case 'T':
		opt->doTestIntra=Yes;
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
	case 'd':
	    opt->noData=true;
	    break;
	case 'K':
		opt->key = optarg;
		break;
	default:
	    cerr << "?? getopt returned character code " <<  c << "\n";
	    break;
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
	   "  --noData            -d  load only one per one fasta file\n"
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

    printSwitch(getSwitch(0));


    getParam(argc,argv,&opt);
    bool error = false;
    bool testResult = true;

    if (opt.doPrintHelp==Yes)
    {
	printHelp();
	exit(0);
    }
    if (opt.doTest || opt.doTestIntra)
    {
    	if (opt.doTest)
        {
    		testResult = executeTests(NB_TEST);


	//		doTest0();
			if(testResult)
			{
				cerr << "Tous les tests sont OK !\n";
				exit(0);
			}
			else
			{
				cerr << "Va vérifier ton code...\n";
				exit(1);
			}
        }
    	if (opt.doTestIntra)
    	{
    		testResult =  testResult && testIntra();
    		testResult = testResult && testIntra2();

    		if(testResult)
			{
				cerr << "Tous les tests sont OK !\n";
				exit(0);
			}
			else
			{
				cerr << "Va vérifier ton code...\n";
				exit(1);
			}
    	}



    }
    else
    {
	if (opt.listFastaPath=="null" && opt.fastaPath=="null")
	{
		if (opt.root=="null")
		{
			if (opt.listFastaPath=="null" && opt.fastaPath=="null")
			{
				cerr << "Need a fasta file (-f file.fasta) or a list of fasta path (-F listFasta_file) \n \tor a path dir (-r pathToTaxon) \n";
				error = true;
			}

		}
		else
		{
			racineBool = true;
		}
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

	FreqKmer *f = NULL;

	if(initFromList)
	{
	    // if(opt.noData)
	    // {
	    // 	cerr << "TOTO\n";
	    // }
	    if(opt.jump==0)
	    {
	    	f = new FreqKmer(opt.windowSize,initFromList,opt.listFastaPath,opt.kmerPath,opt.noData,opt.key);
	    }
	    else
	    {
	    	f = new FreqKmer(opt.windowSize,opt.jump,initFromList,opt.listFastaPath,opt.kmerPath,opt.noData,opt.key);
	    }
	}
	else
	{
	    if(!racineBool)
	    {
	    	if(opt.jump==0)
			{
	    		 f = new FreqKmer(opt.windowSize,initFromList,opt.fastaPath,opt.kmerPath,opt.noData,opt.key);
			}
	    	else
	    	{
	    		f = new FreqKmer(opt.windowSize,opt.jump,initFromList,opt.fastaPath,opt.kmerPath,opt.noData,opt.key);
	    	}
	    }
	    else
	    {
	    	if(opt.jump==0)
			{
	    		f = new FreqKmer(opt.windowSize,opt.kmerPath,opt.noData,opt.root,opt.key);
			}
	    	else
	    	{
	    		f = new FreqKmer(opt.windowSize,opt.jump,opt.kmerPath,opt.noData,opt.root,opt.key);
	    	}
	    }
	}

	f->fillFreq();
	cerr << "Fin fill \n";
	// f->imprimeCSV(opt.output);
	
	delete f;
    }
    return 0;
}
