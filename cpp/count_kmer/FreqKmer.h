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
	Data **data;
	int *indexLineData;
	int **indexLineDataSeq;

	/**
	 * Pour les kmers ## et # le tableau kmerSpace sera
	 * kmerSpace[0]=15; kmerSpace[1]=15+4=19
	 * utile pour enlever un kmer si on enleve le #:
	 * on retire les colonnes de 15+1 à 19.
	 */
	int* kmerSpace;

	int nLigne, nCol, nPattern, tailleFenetre, nData, nbFichierFasta,premier,dernier,index;

public:
	FreqKmer(int tailleF);
	FreqKmer();
	virtual ~FreqKmer();

	void	initPatterns(string fichier);
	void 	initFreq();
	void 	initFromList(string fichier);
	void 	initFromFasta(string fichier);
	void 	copieLigneFreq(int src, int dest,int p);
	void 	compteFenetre2(int *seq,int seq_taille,int debut,int col,Pattern *p);
	void 	compteFenetre(int *seq,int seq_taille,int pos,int indicePattern);
	void 	imprimeCSV(string output);
	void 	add_one(int *seq,int i,int seq_taille,int p);
	void	count(int *seq,int seq_taille,int indicePattern);
	void 	fillFreq();
	int		getCol(int indicePattern,int *seq,int pos);
	double** getFreq(){return freq;}
	int 	getNPattern(){return nPattern;}
	int 	getNbFichierFasta(){return nbFichierFasta;}
	int 	getNCol(){return nCol;}
	int 	getNLigne(){return nLigne;}
	int 	getTailleFenetre(){return tailleFenetre;}
	int 	getNData(){return nData;}
	int 	getStartColKmer(int i);
	int 	getEndColKmer(int i);
	Data** 	getData(){return data;}
	int		getNbLineData(int i);
	int 	getStartLineData(int i);
	int 	getEndLineData(int i);
	int 	getStartLineDataSeq(int i,int j);
	int 	getEndLineDataSeq(int i,int j);
	int 	getNbLineDataSeq(int i,int j);
};

#endif /* FREQKMER_H_ */
