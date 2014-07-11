/*
 * FreqKmer.h
 *
 *  Created on: 10 juil. 2014
 *      Author: jeremy
 */
#include "classPattern.h"
#include "classData.h"

#ifndef FREQKMER_H_
#define FREQKMER_H_

class FreqKmer {

private:
	double **freq;
	Pattern **patterns;
	Data *data;

	int nLigne, nCol, nPattern, tailleFenetre, nData;

public:
	FreqKmer(int tailleF);
	FreqKmer();
	virtual ~FreqKmer();

	void	initPatterns(string fichier);
	void 	initData(string fichier);
	int 	getNPattern(){return nPattern;}
	int 	getNCol(){return nCol;}
	int 	getNLigne(){return nLigne;}
};

#endif /* FREQKMER_H_ */
