/*
 * FreqKmer.cpp
 *
 *  Created on: 10 juil. 2014
 *      Author: jeremy
 */

#include "FreqKmer.h"
#include "classPattern.h"
#include <fstream>
#include <string>

using namespace std;

FreqKmer::FreqKmer() {
	/* On veut les accessions donc Yes */
	data=NULL;
	freq=NULL;
	patterns=NULL;
	nCol=0;
	nData=0;
	nLigne=0;
	nPattern=0;
	tailleFenetre=0;

}

FreqKmer::FreqKmer(int tailleF) {
	data=NULL;
	freq=NULL;
	patterns=NULL;
	nCol=0;
	nData=0;
	nLigne=0;
	nPattern=0;
	tailleFenetre=tailleF;

}

FreqKmer::~FreqKmer() {
	delete data;
	for (int var = 0; var < nPattern ; var++)
	{
		delete patterns[var];
	}
	delete[] patterns;
}

void FreqKmer::initPatterns(string fichier)
{
	int tailleLigne=0;
	string ligne;
	ifstream file(fichier.c_str());
	getline(file,ligne);
	tailleLigne = ligne.length();
    while (file)
    {
    	if(tailleLigne!=0)
    	{
    		nPattern++;
    	}
    	getline(file,ligne);
		tailleLigne = ligne.length();
    }

    ifstream file2(fichier.c_str());

	patterns = new Pattern*[nPattern];
	for (int var = 0; var < nPattern; var++)
	{

		getline(file2,ligne);
		patterns[var] = new Pattern(ligne);

	}

	for(int i=0 ; i<nPattern ; i++)
	{

		nCol += patterns[i]->getAllCombi();

	}

}

void FreqKmer::initFromList(string fichier)
{
	string ligne;
	ifstream file(fichier.c_str());
	getline(file,ligne);
	int nbFichier = 0;
	int tailleLigne = ligne.length();
	while (file)
	{
		if(tailleLigne!=0)
		{
			nbFichier++;
		}
	}
}

void FreqKmer::initFromFasta(string fichier)
{
	data=new Data*[1];
	data[0] = new Data(Yes);
	data[0]->initFrom(fichier,Fasta);
	nData = data[0]->getNtaxa();
//	cout << "ntaxa = " << nData << "\n";
	int tailleSeq=0;
	if (tailleFenetre>0)
	{
		for (int var = 0; var < nData; var++)
		{
//			cout << "FAIL !!! et var = " << var << "\n";
//			cout.flush();
//			cout << " seq = " << (data[0]->getPrimarySequence(var)).length()  << "\n";
			tailleSeq = data[0]->getPrimarySequence(var).length();
//			cout << "FAIL 2 !!! \n";
//			cout.flush();
			if (tailleSeq<tailleFenetre)
			{
				nLigne += 1;
			}
			else
			{
				nLigne += tailleSeq-tailleFenetre+1;
			}
		}
	}
	else
	{
		nLigne = nData;
	}
}


