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

	double **freq; /* tableau des frequences */
	Pattern **patterns; /* les patterns */
	Data **data; /* les sequences */
	int *indexLineData; /* Map pour les Data */
	int **indexLineDataSeq; /* Map pour les seq de Data */

	/**
	 * Pour les kmers ## et # le tableau kmerSpace sera
	 * kmerSpace[0]=15; kmerSpace[1]=15+4=19
	 * utile pour enlever un kmer si on enleve le #:
	 * on retire les colonnes de 15+1 à 19.
	 */
	int *kmerSpace;	/* Map pour les kmers */
	int *tabPremier;
	int *tabDernier;

	int nLigne, nCol, nPattern, tailleFenetre, nData, nbFichierFasta,index,shift;

	/**
	 * Effectue le comptage dans une fenetre
	 * @param	seq: la sequence où effectuer le comptage
	 * @param	seq_taille: taille de la sous sequence (fenetre)
	 * @param	pos: là où commencer le comptage
	 * @param 	indicePattern: indice du pattern courant
	 */
	void 	compteFenetre(int *seq,int seq_taille,int pos,int indicePattern);


	/**
	 * Effectue le comptage de kmer
	 * @param	seq: la sequence où effectuer le comptage
	 * @param	seq_taille: taille de la sous sequence
	 * @param 	indicePattern: indice du pattern courant
	 */
	void	count(int *seq,int seq_taille,int indicePattern);


public:

	/**********************************************************************************/
	/**
	 * Contruit l'object FreqKmer
	 * @param	tailleF: taille de la fenetre
	 * @param	shift: taille du decalage en nucletoide
	 */
	FreqKmer(int tailleF, int shift);

	/**
	 * Contruit l'object FreqKmer
	 * @param	tailleF: taille de la fenetre
	 */
	FreqKmer(int tailleF);

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
	int 		getNLigne(){return nLigne;}
	int 		getTailleFenetre(){return tailleFenetre;}
	int 		getNData(){return nData;}
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
	int			obtainColIndex(int indicePattern,int *seq,int pos);

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
