/*
 * FreqKmer.h
 *
 *  Created on: 10 juil. 2014
 *      Author: jeremy
 */
#include "classPattern.h"
#include "classData.h"
#include "popphyl.h"

#ifndef FREQKMER_H_
#define FREQKMER_H_

class FreqKmer {

private:


	double **freq; /* tableau des frequences de taille nLine*nCol */
	Pattern **patterns; /* tableau de pointeur de pattern définissant les kmers */
	Data **data; /* Tableau de pointeur sur les jeux de données contenant les séquences à analyser */
	int *indexLineData; /* indication sur la ligne de fin dans table freq pour une data donnée */
	int **indexLineDataSeq; /* indication sur la ligne de fin dans table freq pour une sequence d'une data donnée */
	int *kmerSpace;	/* Map pour les kmers avec une taille nPattern, permet de se deplacer horizontalement */


	int nLine, /* Nombre de ligne du tableau freq: nombre de vecteur de frequence */
	nCol, /* Nombre de colonne définit par les patterns: nombre de kmer possible */
	nPattern, /* Nombre de pattern définissant les kmers */
	winSize, /* Taille de fenetre dans laquelle il faut effectuer le comptage de kmer */
	nSeq, /* Nombre total de sequence à traiter */
	nbFichierFasta, /* Nombre de fichier fasta entrée */
	index, /* marqueur de ligne lors du comptage */
	shift; /* taille du decalage lors du passage d'une fenetre à la suivante > 0 */

	/**
	 * Effectue le comptage dans une fenetre
	 * @param	seq: la sequence où effectuer le comptage
	 * @param	win_length: taille de la sous sequence (fenetre)
	 * @param	pos: là où commencer le comptage
	 * @param 	indexPattern: indice du pattern courant
	 */
	void 	winCount(int *seq,int win_length,int pos,int indexPattern,int *current);


	/**
	 * Effectue le comptage de kmer
	 * @param	seq: la sequence où effectuer le comptage
	 * @param	seq_length: taille de la sous sequence
	 * @param 	indexPattern: indice du pattern courant
	 */
	void	count(int *seq,int seq_length,int indexPattern);

	void swapBuffAndCount(int *current,int *previous,int buf_size, int indexPattern,int *seq,int pos);
	void swap(int *current,int *previous,int buf_size);
	void printBuf(int *buf,int buf_size);

public:

	/**********************************************************************************/
	/**
	 * Contruit l'object FreqKmer
	 * @param	win_size: taille de la fenetre
	 * @param	shift: taille du decalage en nucletoide
	 */
	FreqKmer(int win_size, int shift);

	/**
	 * Contruit l'object FreqKmer
	 * @param	win_size: taille de la fenetre
	 */
	FreqKmer(int win_size);

	FreqKmer();
	virtual ~FreqKmer();

	/**********************************************************************************/

	/**
	 * Initialise la table patterns
	 * @param	listPatterns fichier contenant les patterns
	 */
	void	initPatterns(string listPatterns);

	/**
	 * Met chaque case
	 * du tableau freq à 0.
	 */
	void 	initFreq();

	/**
	 * Initialise la table Data
	 * @param 	listFasta fichier contenant les chemins des fichiers fasta à init
	 */
	void 	initDataFromListFastaPath(string listFasta);

	/**
	 * Initialise la table Data
	 * @param 	fasta fichier contenant les sequences
	 */
	void 	initFromFasta(string fasta);

	/**
	 * Imprime au format csv la table freq
	 * @param 	output fichier où ecrire le csv
	 */
	void 	imprimeCSV(string output);

	/**
	 * Methode principale qui permet
	 * de remplir le tableau
	 * après avoir chargé les données
	 */
	void 		fillFreq();

	/*******************************ACCESSEUR*****************************************/
	/**
	 * inline accesor
	 */
	double** 	getFreq(){return freq;}
	int 		getNPattern(){return nPattern;}
	int 		getNbFichierFasta(){return nbFichierFasta;}
	int 		getNCol(){return nCol;}
	int 		getNLine(){return nLine;}
	int 		getWinSize(){return winSize;}
	int 		getNSeq(){return nSeq;}
	Data** 		getData(){return data;}


	/*******************************DEPLACEMENT HORIZONTAL*****************************************/
	/**
	 * @param i: indice du kmer
	 * @return l'indice de début de colonne
	 * du kmer i
	 */
	int 		obtainStartColKmer(int i);

	/**
	 * @param i: indice du kmer
	 * @return l'indice de fin de colonne
	 * du kmer i
	 */
	int 		obtainEndColKmer(int i);

	/**
	 * @param indicePatter: indice du kmer
	 * @param seq: la sequence courante
	 * @param pos: index dans la sequence ou recuperer le kmer
	 * @return l'indice de fin de colonne
	 * du kmer i
	 */
	int			obtainColIndex(int indexPattern,int *seq,int pos);

	/*******************************DEPLACEMENT VERTICAL*****************************************/
	/******** Pour un Data ********/
	/**
	 * @param i: indice de l'objet Data
	 * @return le nombre de ligne qu'occupe Data[i]
	 */
	int			obtainNbLineData(int i);

	/**
	 * @param i: indice de l'objet Data
	 * @return l'indice de la ligne où commence
	 * les fréquences pour Data[i]
	 */
	int 		obtainStartLineData(int i);

	/**
	 * @param i: indice de l'objet Data
	 * @return l'indice de la ligne où termine
	 * les fréquences pour Data[i]
	 */
	int 		obtainEndLineData(int i);

	/******** Pour une sequence d'un data ********/

	/**
	 * @param i: indice de l'objet Data
	 * @param j: indice de la sequence dans Data[j]
	 * @return l'indice de la ligne où commence
	 * les frequences de la j-ème
	 * sequence de Data[i]
	 */
	int 		obtainStartLineDataSeq(int i,int j);

	/**
	 * @param i: indice de l'objet Data
	 * @param j: indice de la sequence dans Data[j]
	 * @return l'indice de la ligne où termine
	 * les frequences de la j-ème
	 * sequence de Data[i]
	 */
	int 		obtainEndLineDataSeq(int i,int j);

	/**
	 * @param i: indice de l'objet Data
	 * @param j: indice de la sequence dans Data[j]
	 * @return le nombre de ligne qu'occupe la j-ème
	 * sequence de Data[i]
	 */
	int 		obtainNbLineDataSeq(int i,int j);

	/**
	 * Permet d'avoir le nombre de decalage entre deux position
	 * et une fenetre et un pas
	 * @param	i: position de depart
	 * @param	j: position d'arrivée
	 * @param	l: taille de la fenetre
	 * @param	z: le pas de decalage
	 * @return le nombre de decalage possible entre i, j avec
	 * une taille de fenetre l et un decalage de z
	 */
	int 		obtainNbLineWindow(int i,int j,int l, int z);

};

#endif /* FREQKMER_H_ */
