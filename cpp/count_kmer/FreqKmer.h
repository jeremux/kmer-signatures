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

/*************************************************************************************************************************
 * *************************************************************************************************************************
 * *****************************************       PRIVATE      ************************************************************
 * *************************************************************************************************************************
 * *************************************************************************************************************************
 */
private:


	double **freq; /* tableau des frequences de taille nLine*nCol */
	Pattern **patterns; /* tableau de pointeur de pattern définissant les kmers */
	Data **data; /* Tableau de pointeur sur les jeux de données contenant les séquences à analyser */
	int *indexLineData; /* indication sur la ligne de fin dans table freq pour une Data donnée */
	int **indexLineDataSeq; /* indication sur la ligne de fin dans table freq pour une sequence d'une data donnée */
	int *kmerSpace;	/* Map pour les kmers avec une taille nPattern, permet de se deplacer horizontalement */
	bool **mask; /* mask pour savoir quels séquences de quel jeu de donnée on considère, de taille nData*nSeq */
	string **taxaDataSeq;
	int *indexTaxaInFasta;

	int nLine, /* Nombre de ligne du tableau freq: nombre de vecteur de frequence */
	nCol, /* Nombre de colonne définit par les patterns: nombre de kmer possible */
	nPattern, /* Nombre de pattern définissant les kmers */
	winSize, /* Taille de fenetre dans laquelle il faut effectuer le comptage de kmer */
	nSeq, /* Nombre total de sequence à traiter */
	nbFastaFile, /* Nombre de fichier fasta entrée */
	shift, /* taille du decalage lors du passage d'une fenetre à la suivante > 0 */
	nbChildTaxa; /* nombre de taxon fils */

	bool noData; /* boolean pour savoir si on instancie les objets data d'une traite */
	string listData; /* fichier de la list des chemins fasta */
	string *pathFasta;
	string pathRoot; /* dossier ou on devra prendre une décision */
	vector<string> pathChildTaxa; /* chemin des taxon fils du dossier courant */
	vector<string> idTaxa; /* id des taxon fils du dossier courant */
	vector<string> listPathFasta;
	/**
	 * Effectue le comptage dans une fenetre
	 * @param	seq: la sequence où effectuer le comptage
	 * @param	win_length: taille de la sous sequence (fenetre)
	 * @param	pos: là où commencer le comptage
	 * @param 	indexPattern: indice du pattern courant
	 * @param	current: buffer qui va garder les frequences trouvés pour la fenetre courante
	 * @param	start_line: là ligne où commence l'écriture dans la table freq
	 */
	void 	winCount(int *seq,int win_length,int pos,int indexPattern,int *current,int start_line);


	/**
	 * Effectue le comptage de kmer
	 * @param	seq: la sequence où effectuer le comptage
	 * @param	seq_length: taille de la sequence
	 * @param 	indexPattern: indice du pattern courant
	 * @param	start_line: là ligne où commence l'écriture dans la table freq
	 */
	void	count(int *seq,int seq_length,int indexPattern,int start_line);

	/**
	 * Recupère la partie commune de previous pour le comptage
	 * et compte les nouveaux kmers rencontrés lors du décalage de la fenetre
	 * @param 	current: la ligne de comptage actuelle
	 * @param	previous: la ligne de comptage effectué avant le decalage (null si premier comptage)
	 * @param 	buf_size: la taille des buffers current et previous
	 * @param	pos: position où on se trouve dans la sequence
	 * @param	seq: la sequence où effectuer le comptage
	 * @param 	indexPattern: indice du pattern courant
	 */
	void copyBuffAndCount(int *current,int *previous,int buf_size, int indexPattern,int *seq,int pos);

	/**
	 * Permet d'intervertir deux buffers
	 * @param	current
	 * @param	previous
	 * @param	buf_size : taille des buffers
	 */
	void swap(int *current,int *previous,int buf_size);

	/**
	 * Permet d'imprimer un buffer sur stdout : [x][y]...
	 * Utilisé essentiellement pour vérifier le bon fonctionnement
	 * @param	buf: buffer à afficher
	 * @param	buf_size: taille du buffer
	 */
	void printBuf(int *buf,int buf_size);

	void initPathFasta(string file);

	void initTabIndexTaxaInFasta(string key_fasta);


/*************************************************************************************************************************
 * *************************************************************************************************************************
 * *****************************************       PUBLIC      ************************************************************
 * *************************************************************************************************************************
 * *************************************************************************************************************************
 */
public:

	/**********************************************************************************/

	/**
	 * Contruit l'object FreqKmer à partir d'un fichier, contenant
	 * la liste des chemins des fichiers fasta à analyser si bool = true
	 * ou à aprtir d'un fichier fasta sinon.
	 * @param	win_size: taille de la fenetre
	 * @param	list: true si le fichier en parametre est une liste de path fasta
	 * @param	file: fichier fasta ou fichier contenant la liste de fichier fasta
	 * @param 	patternFile: fichier contenant les patterns de kmer
	 * @param	noData	: true si on ne charge pas toutes les données d'un coup
	 */
	FreqKmer(int win_size,bool list, string file,string patternFile, bool noData,string pathRoot);

	/**
	 * Contruit l'object FreqKmer à partir d'un fichier, contenant
	 * la liste des chemins des fichiers fasta à analyser
	 * @param	win_size: taille de la fenetre
	 * @param	shift: taille du decalage en nucletoide
	 * @param	fastaFile: fichier fasta
	 * @param 	patternFile: fichier contenant les patterns de kmer
	 * @param	noData	: true si on ne charge pas toutes les données d'un coup
	 */
	FreqKmer(int win_size,int shift,bool list, string file,string patternFile, bool noData,string pathRoot);


	FreqKmer();

	virtual ~FreqKmer();

	/**********************************************************************************/
	bool	takeDataSeq(int indexData,int indexSeq);
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
	 * Initialise la table Data. Lit les fastas définit
	 * par leur chemin dans listFasta. Met à jour le nombre de ligne nLigne
	 * et le nombre de sequence nSeq et le nombre de fichier fasta nbFichierFasta
	 * @param 	listFasta fichier contenant les chemins des fichiers fasta à init
	 */
	void 	initDataFromListFastaPath(string listFasta);

	/**
	 * Initialise la table Data. Lit les fastas définit
	 * par leur chemin dans listFasta. Met à jour le nombre de ligne nLigne,
	 * le nombre de sequence nSeq et le nombre de fichier fasta nbFichierFasta (= 1)
	 * @param 	fasta fichier fasta contenant les données
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
	 * Pour la première fenetre de la
	 * première sequence du premier jeu de donnée
	 * Va ecrire à la premiere ligne
	 * 		à la premiere colonne: le nombre du premier kmer trouvé dans la fenetre
	 * 		à la seconde  colonne: le nombre du premier kmer trouvé dans la fenetre
	 * 		...
	 *
	 * ....
	 * Pour la derniere fenetre de la
	 * derniere sequence du dernier jeu de donnée
	 * Va ecrire à la derniere ligne
	 * 		à la premiere colonne: le nombre du premier kmer trouvé dans la fenetre
	 * 		à la seconde  colonne: le nombre du premier kmer trouvé dans la fenetre
	 * 		...
	 */
	void 		fillFreq();

	int 		getDirTaxonFromPath(string path,int &nbChild, vector<string> &files,vector<string> &taxids);

	/*******************************ACCESSEUR*****************************************/
	/**
	 * inline accesor
	 */
	double** 	getFreq(){return freq;}
	int 		getNPattern(){return nPattern;}
	int 		getNbFichierFasta(){return nbFastaFile;}
	int 		getNCol(){return nCol;}
	int 		getNLine(){return nLine;}
	int 		getWinSize(){return winSize;}
	int 		getNSeq(){return nSeq;}
	int			getNbChild(){return nbChildTaxa;}
	Data** 		getData(){return data;}


	/*******************************DEPLACEMENT HORIZONTAL*****************************************/
	/**
	 * Permet de récupérer le numéro de colonne
	 * où commence un kmer dans la table freq
	 * @param i: indice du kmer
	 * @return l'indice de début de colonne
	 * du kmer i
	 */
	int 		obtainStartColKmer(int i);

	/**
	 * Permet de récupérer le numéro de colonne
	 * où termine un kmer dans la table freq
	 * @param i: indice du kmer
	 * @return l'indice de fin de colonne
	 * du kmer i
	 */
	int 		obtainEndColKmer(int i);

	/**
	 * Permet de savoir quel kmer incrémenter
	 * @param indexPattern: indice du kmer
	 * @param seq: la sequence courante
	 * @param pos: index dans la sequence ou recuperer le kmer
	 * @return la position du kmer dans la table freq
	 */
	int			obtainColIndex(int indexPattern,int *seq,int pos);



	/*******************************DEPLACEMENT VERTICAL*****************************************/
	/******** Pour un Data ********/
	/**
	 * Permet de savoir le nombre de ligne
	 * utilisé par un jeu de données
	 * @param i: indice de l'objet Data
	 * @return le nombre de ligne qu'occupe Data[i]
	 */
	int			obtainNbLineData(int i);

	/**
	 * Peremet de savoir où commence
	 * les fréquences d'un jeu de données dans
	 * la table freq
	 * @param i: indice de l'objet Data
	 * @return l'indice de la ligne où commence
	 * les fréquences pour Data[i]
	 */
	int 		obtainStartLineData(int i);

	/**
	 * Peremet de savoir où termine
	 * les fréquences d'un jeu de données dans
	 * la table freq
	 * @param i: indice de l'objet Data
	 * @return l'indice de la ligne où termine
	 * les fréquences pour Data[i]
	 */
	int 		obtainEndLineData(int i);

	/******** Pour une sequence d'un data ********/

	/**
	 * Peremet de savoir où commence
	 * les fréquences d'une sequence d'un jeu de données
	 * précis dans la table freq
	 * @param i: indice de l'objet Data
	 * @param j: indice de la sequence dans Data[j]
	 * @return l'indice de la ligne où commence
	 * les frequences de la j-ème
	 * sequence de Data[i]
	 */
	int 		obtainStartLineDataSeq(int i,int j);


	/**
	 * Peremet de savoir où termine
	 * les fréquences d'une sequence d'un jeu de données
	 * précis dans la table freq
	 * @param i: indice de l'objet Data
	 * @param j: indice de la sequence dans Data[j]
	 * @return l'indice de la ligne où termine
	 * les frequences de la j-ème
	 * sequence de Data[i]
	 */
	int 		obtainEndLineDataSeq(int i,int j);



	/**
	 * Peremet de savoir le nombre de ligne
	 * des fréquences d'une sequence d'un jeu de données
	 * précis dans la table freq
	 * @param i: indice de l'objet Data
	 * @param j: indice de la sequence dans Data[j]
	 * @return le nombre de ligne qu'occupe la j-ème
	 * sequence de Data[i]
	 */
	int 		obtainNbLineDataSeq(int i,int j);

	int 		obtainStartLineTaxaInFastaList(int i);
	int 		obtainEndLineTaxaInFastaList(int i);
	int 		obtainNbLineTaxaInFastaList(int i);

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

	void		obtainDataSeqFromLine(int line,int &idData,int &idSeq);

	void 		obtainDataSeqWinFromLine(int line,int &idData,int &idSeq,int &idWin);

	string 		getPathChildTaxa(int i);

	string 		getIdTaxa(int i);

	void 		writeListFasta();

};

#endif /* FREQKMER_H_ */
