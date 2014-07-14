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


#include "classPattern.h"
#include "classData.h"
#include "FreqKmer.h"


int main(int argc, char **argv) {

//	Data *d = new Data(Yes);
//	d->initFrom("file.fasta",Fasta);


	FreqKmer *f = new FreqKmer(4);
	string filename = "liste.txt";
//	cout << filename << "\n";
	f->initPatterns("pattern.txt");
	Pattern *p = new Pattern("##");


//	f->initFromFasta("file2.fasta");
	f->initFromList(filename);


	//int nbTaxa = d->getNtaxa();
//	cout << "sequence = \n" << seq <<"\n";
//	cout << "accession = " << acc << "\n";
//	cout << "Nombre de sequences = " << nbTaxa << "\n";

	f->fillFreq();
	f->imprimeCSV();
	return 0;
}





