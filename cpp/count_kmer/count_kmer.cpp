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


	FreqKmer *f = new FreqKmer(50);
	string filename = "/home/jeremy/mitomer/trunk/create_db/Eukaryota__2759/";
	filename += argv[1];
	cout << filename << "\n";
	f->initPatterns("pattern.txt");


//	f->initFromFasta("file2.fasta");
	f->initFromFasta(filename);


	//int nbTaxa = d->getNtaxa();
//	cout << "sequence = \n" << seq <<"\n";
//	cout << "accession = " << acc << "\n";
//	cout << "Nombre de sequences = " << nbTaxa << "\n";

//	int nbPattern = f->getNPattern();
//	int nbCol = f->getNCol();
//	cout << "Nombre de Pattern = " << nbPattern << "\n";
//	cout << "Nombre de col = " << nbCol << "\n";
	cout << f->getNLigne() << "\n";
	return 0;
}





