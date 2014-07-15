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
//	cout << "Debut initPatterns\n";
	f->initPatterns("pattern.txt");
//	cout << "Fin initPatterns\n";

//	f->initFromFasta("file2.fasta");
//	cout << "Debut init fasta\n";
	f->initFromList(filename);
//	cout << "Fin init fasta\n";


	//int nbTaxa = d->getNtaxa();
//	cout << "sequence = \n" << seq <<"\n";
//	cout << "accession = " << acc << "\n";
//	cout << "Nombre de sequences = " << nbTaxa << "\n";
//	int col = f->getCol(1,f->getData()[0]->getDataObject()[0],2);
	//cerr << "Col recherchee = " << col << "\n";
//	cout << "Debut fillFreq\n";
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
//	freq[0][16] = 2;
//	freq[0][17] = 2;
//	freq[1][16] = 1;
//	freq[1][17] = 2;
//	freq[1][19] = 1;
//	freq[2][16] = 1;
//	freq[2][17] = 1;
//	freq[2][18] = 1;
//	freq[2][19] = 1;
//	freq[3][16] = 1;
//	freq[3][17] = 1;
//	freq[3][19] = 2;
//	freq[4][16] = 1;
//	freq[4][17] = 1;
//	freq[4][19] = 2;
//	freq[5][16] = 1;
//	freq[5][17] = 1;
//	freq[5][18] = 1;
//	freq[5][19] = 1;
//	freq[6][16] = 1;
//	freq[6][17] = 2;
//	freq[6][18] = 1;

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
		cerr << "Success !\n";
	}
	else
	{
		cerr << "Try again...\n";
	}

//	for(int j=0;j<col;j++)
//		{
//			cout << j << ";";
//		}
//		cout << "\n";
//		for(int i=0;i<ligne;i++)
//		{
//			for(int j=0;j<col;j++)
//			{
//				cout << freq[i][j] << ";";
//			}
//			cout << "\n";
//		}
//	f->imprimeCSV();
//	cout << "Fin Impression\n";
	return 0;
}





