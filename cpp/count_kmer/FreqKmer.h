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

/**
 * @class FreqKmer
 * @brief Permet d'evaluer la frequence en k-mer
 */
class FreqKmer {

/*************************************************************************************************************************
 * *************************************************************************************************************************
 * *****************************************       PRIVATE      ************************************************************
 * *************************************************************************************************************************
 * *************************************************************************************************************************
 */
private:


	double **freq; /**<  tableau des frequences de taille nLine*nCol */
	Pattern **patterns; /**<  tableau de pointeur de pattern définissant les kmers */
	Data **data; /**<  Tableau de pointeur sur les jeux de données contenant les séquences à analyser */
	int *indexLineData; /**<  indication sur la ligne de fin dans table freq pour une Data donnée */
	int **indexLineDataSeq; /**<  indication sur la ligne de fin dans table freq pour une sequence d'une data donnée */
	int *kmerSpace;	/**<  Map pour les kmers avec une taille nPattern, permet de se deplacer horizontalement */
	bool **mask; /**<  mask pour savoir quels séquences de quel jeu de donnée on considère, de taille nData*nSeq */
	string **taxaDataSeq;/**< Nom du taxon pour une sequence d'un taxon */
	int *indexTaxaInFasta; /**< Indice de fin du taxa dans le fichier fasta */
	bool freqFilled; /**< true si la table freq est entièrement remplie */
	bool isList; /**< Si le fichier d'entrée est une liste de fasta */
	string pathFastaFile;/**< Chemin du fichier fasta */
	string pathPattern; /**< Chemin du fichier contenant la liste des patterns de kmer @see classPattern.h */
	string key_fasta; /**< cox1, cox2, cox3, cytb ou genomes */
	bool initFromRoot; /**< Initialisation depuis un dossier taxa */
	bool initWithJump; /**< Le decalage pour le comptage est fourni */
	int *nbSeq; /**< Le nombre de sequence par jeu de données */


	int nLine, /**<  Nombre de ligne du tableau freq: nombre de vecteur de frequence */
	nCol, /**<  Nombre de colonne définit par les patterns: nombre de kmer possible */
	nPattern, /**<  Nombre de pattern définissant les kmers */
	winSize, /**<  Taille de fenetre dans laquelle il faut effectuer le comptage de kmer */
	nSeq, /**<  Nombre total de sequence à traiter */
	nbFastaFile, /**<  Nombre de fichier fasta entrée */
	shift, /**<  taille du decalage lors du passage d'une fenetre à la suivante > 0 */
	nbChildTaxa; /**<  nombre de taxon fils */

	bool noData; /**<  boolean pour savoir si on instancie les objets data d'une traite */
	string listData; /**<  fichier de la list des chemins fasta */
	string *pathFasta;/**< La liste des chemins vers les fichiers fasta */
	string pathRoot; /**<  dossier ou on devra prendre une décision */
	vector<string> pathChildTaxa; /**<  chemin des taxon fils du dossier courant */
	vector<string> idTaxa; /**<  id des taxon fils du dossier courant */
	vector<string> listPathFasta;/**< La liste des chemins vers les fichiers fasta pour un taxa */
	vector<string> idTaxaFromData; /**< Liste des taxid dans l'ordre de @see listPathFasta */
	vector<pair<int, int> > sampledTaxon; /**< Les sequences après echantillonage */

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


	/**
	 * Initialise les chemins fasta
	 * lorsqu'une liste est fourni
	 * @param file: chemin du fichier fasta
	 */
	void initPathFasta(string file);


	/**
	 * Initialise la liste de fichier fasta
	 * des sous taxons du taxon (dossier) courant
	 * @param key_fasta:	mot-clé de la sequence a recuperer
	 * 						cox1,cox2,...,genomes.
	 */
	void initTabIndexTaxaInFasta(string key_fasta);

	/**
	 * Methode permettant de "cacher" une sequence
	 * @param i:	indice du jeu de donnee a cacher
	 * @param j:	indice de la sequence du i-ème jeu de donnes
	 */
	void setFalseMask(int i,int j);
	/**
	 * Permet de normaliser une
	 * ligne du tableau freq
	 * @param indexPattern: indice du pattern
	 * @param lin
	 */
	void normalizeLine(int indexPattern,int line);

	/**
	 * Ecrit l'en tête du fichier
	 * au foramt weka
	 * @param os: ofstream d'écriture
	 */
	void writeHeaderWeka(ofstream &os);

	/**
	 * Ecrit la ligne de frequence
	 * au format weka
	 * @param os: ofstream d'écriture
	 * @param i: indice du jeu de données
	 * @param j: j-ème sequence de data[j]
	 */
	void writeLineInOs(ofstream &os,int i,int j);

	/**
	 * Permet de recuperer un taxid
	 * avec la struture de path definit
	 * par la construction de la base de
	 * données locale : taxon_taxid1/taxon_taxid2/...
	 * @param line: chemin du taxon
	 * @return le taxid du chemin line
	 */
	string getTaxidFromString(string line);

	/**
	 * Ecrit la validation croisée
	 * au foramt weka
	 * @param percent: pourcentage pour la prediction
	 * @param outLearn: fichier de sortie pour l'apprentissage
	 * @param outPredict: fichier de sortie pour la prediction
	 * @see writeLineInOs
	 */
	void writeCrossVal(int percent,string outLearn,string outPredict);

	/**
	 * Ecrit la ligne de frequence
	 * au format weka
	 * @param os: ofstream d'écriture
	 * @param i: indice du jeu de données
	 * @param j: j-ème sequence de data[j]
	 * @param f: objet à consiférer
	 */
	void writeLineInOs(ofstream &os,int i,int j,FreqKmer *f);

	/**
	 * Méthode permetttant
	 * de changer le mode noData en
	 * cours de route
	 * @param b: nouvelle valeure de noData
	 */
	void setNoData(bool b);

	/**
	 * Méthode permettant
	 * de charger toutes les données
	 * notamment après un
	 * setNoData(true)
	 */
	void loadAll();

	/**
	 * méthode permettant de savoir
	 * si dans un jeu de données
	 * on considère au moins une séquences
	 * @param i: indice du jeu de données à tester
	 * @return : true si au moins une séquences est à
	 * 			 considérer dans data[i]
	 */
	bool existSeqTrueInData(int i);


/*************************************************************************************************************************
 * *************************************************************************************************************************
 * *****************************************       PUBLIC      ************************************************************
 * *************************************************************************************************************************
 * *************************************************************************************************************************
 */
public:

	/**********************************************************************************/

	/**
	 * Contruscteur le decalage = 20% de win_size
	 * @param win_size:   taille de la fenetre
	 * @param list:       booleen, false si le fichier file en parametre et un fichier fasta
	 * 					           true  si le fichier file en parametre est une liste de fichier fasta
	 * @param file:       chemin du fichier fasta ou du fichier de la liste des chemins fasta
	 * @param patternFile chemin du fichier contenant les patterns de kmer
	 * @param b:          noData, si true alors les données ne sont chargees que lors de leur utilisation
	 * @param key:        mot cles pour le comptage de frequence (cox1,cox2,...,genomes)
	 */
	FreqKmer(int win_size,bool list, string file,string patternFile, bool noData,string key);

	/**
	 * Contruscteur
	 * @param win_size:   taille de la fenetre
	 * @param s:          taille du decalage
	 * @param list:       booleen, false si le fichier file en parametre et un fichier fasta
	 * 					           true  si le fichier file en parametre est une liste de fichier fasta
	 * @param file:       chemin du fichier fasta ou du fichier de la liste des chemins fasta
	 * @param patternFile chemin du fichier contenant les patterns de kmer
	 * @param b:          noData, si true alors les données ne sont chargees que lors de leur utilisation
	 * @param key:        mot cles pour le comptage de frequence (cox1,cox2,...,genomes)
	 */
	FreqKmer(int win_size,int s,bool list, string file,string patternFile, bool b,string key);

	/**
	 * Contruscteur le decalage = 20% de win_size
	 * @param win_size:   taille de la fenetre
	 * @param patternFile chemin du fichier contenant les patterns de kmer
	 * @param b:          noData, si true alors les données ne sont chargees que lors de leur utilisation
	 * @param pathR: 	  chemin du taxon où établir l'apprentissage
	 * @param key:        mot cles pour le comptage de frequence (cox1,cox2,...,genomes)
	 */
	FreqKmer(int win_size,string patternFile, bool b,string pathR,string key);

	/**
	 * Contruscteur
	 * @param win_size:   taille de la fenetre
	 * @param s:		  taille du decalage
	 * @param patternFile chemin du fichier contenant les patterns de kmer
	 * @param b:          noData, si true alors les données ne sont chargees que lors de leur utilisation
	 * @param pathR: 	  chemin du taxon où établir l'apprentissage
	 * @param key:        mot cles pour le comptage de frequence (cox1,cox2,...,genomes)
	 */
	FreqKmer(int win_size,int s,string patternFile, bool b,string pathR,string key);



	/**
	 *
	 */
	virtual ~FreqKmer();

	/**
	 * Permet d'obtenier un echantillon du jeu
	 * de donnees initialement fourni
	 * @param list:	liste de couple d'entier des séquences
	 * 				à considérer
	 * return:				un pointeur d'objet de type FreqKmer
	 */
	FreqKmer* sampleMe(vector<pair<int, int> > list);



	typedef std::pair<int, int> intPair;
	/**********************************************************************************/

	/**
	 * Methode pour savoir si on
	 * considere une sequence dans nos calcul
	 * @param i:	indice du jeu de donnee a questionner
	 * @param j:	indice de la sequence du i-ème jeu de donnes
	 */
	bool	takeDataSeq(int indexData,int indexSeq);

	/**
	 * Permet d'initialiser le tableau
	 * patterns[] a partir du fichier en parametre
	 * @param fichier	chemin du fichier de patterns de kmer
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

	/**
	 * Methode permettant de recuperer la liste
	 * des sous taxons du taxon courant
	 * @param dir:		chemin du taxon courant
	 * @param nbChild:	où sauvegarder le nombre de sous taxon
	 * @param files:	liste de string ou sera enregistré les chemins vers les fastas
	 * @param taxids:	liste de string contenant les taxids des sous taxons du taxon courant
	 */
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
	string 		getpathFastaFile(){return pathFastaFile;}
	vector<pair<int, int> > getSampledTaxon(){return sampledTaxon;}


	/*******************************DEPLACEMENT HORIZONTAL*****************************************/
	/**
	 * Permet de se deplacer horizontalement
	 * @param i:	index du pattern
	 * return:	    l'indice du début de la colonne du i-ème pattern
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


	/**
	 * Permet d'avoir le numéro de ligne (qui
	 * commence à 0) de début
	 * dans la liste des chemin fasta pour un taxon donnée
	 * @param i:	index du taxon
	 * return:		l'indice de début de ligne du taxon i
	 */
	int 		obtainStartLineTaxaInFastaList(int i);

	/**
	 * Permet d'avoir le numéro de ligne (qui
	 * commence à 0) de fin
	 * dans la liste des chemin fasta pour un taxon donnée
	 * @param i:	index du taxon
	 * return:		l'indice de fin de ligne du taxon i
	 */
	int 		obtainEndLineTaxaInFastaList(int i);

	/**
	 * Permet d'avoir le nombre de fichier fasta
	 * pour un taxon donnée
	 * @param i:	index du taxon
	 * return:		le nombre de fichier fasta pour le i-ème taxon
	 */
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

	/**
	 * Methode pour recuperer l'indice de la sequence et du jeu de
	 * donnee à partir du ligne de la table freq
	 * @param line:		indice de ligne dans la table freq
	 * @param idData:	où sauvegarder l'indice du jeu de donnee
	 * @param idSeq:	où sauvegarder l'indice de la sequence
	 */
	void		obtainDataSeqFromLine(int line,int &idData,int &idSeq);

	/**
	 * Methode pour recuperer l'indice de la sequence et du jeu de
	 * donnee à partir du ligne de la table freq et egalement la fenetre
	 * @param line:		indice de ligne dans la table freq
	 * @param idData:	où sauvegarder l'indice du jeu de donnee
	 * @param idSeq:	où sauvegarder l'indice de la sequence
	 * @param idWin: 	où sauvegarder l'indice de la fenetre
	 */
	void 		obtainDataSeqWinFromLine(int line,int &idData,int &idSeq,int &idWin);


	/**
	 * Permet de recupere le chemin du
	 * dossier d'indice dir_i
	 * @param dir_i:	indice du dossier
	 * return: 			le chemin du dossier d'indice dir_i
	 */
	string 		getPathChildTaxa(int i);

	/**
	 * Permet de recuperer le taxid
	 * à partir du dossier d'indice dir_i
	 * @param dir_i:	indice du dossier
	 * return:			le taxid du dossier d'indice dir_i
	 */
	string 		getIdTaxa(int i);


	/**
	 * Methode permettant d'ecrire en
	 * dur dans le dossier courant la liste des fasta
	 * des sous taxon du dossier courant
	 */
	void 		writeListFasta();

	/**
	 * Tirage aleatoire sans remise
	 * @param result:		liste des entiers tirés
	 * @param tabSize:		borne sup du tirage
	 * @param sampleSize:	nombre d'entier à tirer
	 */
	void randomTab(vector<int> *result,int tabSize,int sampleSize);

	/**
	 * Permet d'obtenier un echantillon du jeu
	 * de donnees initialement fourni
	 * @param sampleSize:	taille de l'échantillon
	 * return:				un pointeur d'objet de type FreqKmer
	 */
	FreqKmer* sampleMe(int sampleSize);

	/**
	 * Methode permettant d'avoir le
	 * nombre de case à true dans un tableau de booleens entre
	 * deux indices
	 * @param tab:		tableu de booleens
	 * @param start:	indice de debut ou compter
	 * @param end:		indice de fin ou arreter le comptage
	 * return 			le nombre de case a true entre les indices
	 * 					start et end incluses dans le tableau tab
	 */
	int 	  getNbTrue(bool *tab,int start,int end);

	/**
	 * Methode permettant de savoir combien
	 * de sequences sont considérées, n'est
	 * pas utile au programme mais permet des tests
	 * basiques
	 */
	int 	  getNbAllTrue();

	/**
	 * @Biref Methode permettant d'obtenir
	 * le nombre de sequence par taxon
	 * @param i: index du taxon
	 * @return : le nombre de sequences dans le taxon i.
	 */
	int		  getNSeqInTaxa(int i);

	/**
	 * @Biref methode permettant de changer les
	 * sequences à considerer
	 * @param candidate: vecteur contenant une lsite d'indice
	 * 				     correspondant aux sequences à garer
	 * @param mask_tmp: tableau a masquer
	 * @param indexTaxa: indice du taxon dont les sequences sont a masquer
	 */
	void 	  maskTab(vector<int> *candidate,bool **mask_tmp, int indexTaxa);

	/**
	 * @Brief methode permettant d'avoir
	 * la somme d'une ligne dans la table @see freq
	 * @param indexPattern: indice du pattern
	 * @param line: indice de la ligne à calculer
	 *
	 */
	double getSum(int indexPattern,int line);

	/**
	 * @Biref normalise les valeurs dans
	 * la table freq
	 */
	void 	  normalize();

	/**
	 * @Brief permet d'ecrire le fichier
	 * de configuration de l'objet courant
	 * @param output: nom du fichier de sortie
	 * @see initFromConf
	 */
	void 	  writeConfFeq(string output);

	/**
	 * Permet de s'initialiser à partir
	 * d'un fichier de configuration
	 * @param conf: fichier de configuration
	 */
	FreqKmer* 	  initFromConf(string conf);

	/**
	 * Teste l'égalité entre
	 * les tableaux freq
	 * @param f: objet FreqKmer a tester
	 * @return : vrai si les tableaux de frequences de f et
	 * 			 et de l'objet courant sont identiques
	 */
	bool 	  equal(FreqKmer *f);

	/**
	 *	Etablit une validation croisée
	 *	pour l'objet courant, ecrit
	 *	en sortie deux fichier learn.arff et
	 *	toPredict.arff
	 *	@param percent: pourcentage pour la prédiction
	 */
	void 	  writeCrossVal(int percent);

	/**
	 * Calcul la fréquence en kmer pour
	 * une sequence
	 * @param data_i: indice du data
	 * @param seq_j: indice de la sequence
	 */
	void	fillFreq(int data_i,int seq_j);

	/**
	 *	Etablit une validation croisée
	 *	pour l'objet courant, ecrit
	 *	en sortie deux fichier learn-i.arff et
	 *	toPredict-i.arff
	 *	@param percent: pourcentage pour la prédiction
	 *	@param i: indice a ecrire dans le nom de fichier
	 */
	void 	writeCrossVal(int percent,int i);

	/**
	 * Etablit une validation croisée entre deux objet
	 * passé en paramètre
	 * @param freqLearn: objet pour l'apprentissage
	 * @param freqPredict: objet pour la prédiction
	 * @param percent: pourcentage pour la prédiction
	 * @param outLearn: nom du fichier de sortie à apprendre
	 * @param outToclassify: nom du fichier de sortie à prédire
	 */
	void 	writeCrossVal(FreqKmer *freqLearn, FreqKmer *freqPredict, int percent, string outLearn, string outToclassify);

	/**
	 * Etablit plusieurs validations croisées entre deux objet
	 * passé en paramètre
	 * @param freqLearn: objet pour l'apprentissage
	 * @param freqPredict: objet pour la prédiction
	 * @param percent: pourcentage pour la prédiction
	 * @param i: nombre de validation croisée à produire
	 */
	void 	writeNCrossVal(FreqKmer *freqLearn, FreqKmer *freqPredict, int percent, int i);

	/**
	 * Etablit plusieurs validations croisées entre deux objet
	 * passé en paramètre
	 * @param freqLearn: objet pour l'apprentissage
	 * @param freqPredict: objet pour la prédiction
	 * @param percent: pourcentage pour la prédiction
	 * @param i: nombre de validation croisée à produire
	 * @param id: id pour suffixer les fichiers de sortie
	 * 			  pour la prédiction.
	 */
	void   writeNCrossVal(FreqKmer *freqLearn, FreqKmer *freqPredict, int percent, int i,string id);

	/**
	 * Etablit une validation croisée entre l'objet courant
	 *  et l'objet passé en paramètre
	 * @param freqPredict: objet pour la prédiction
	 * @param percent: pourcentage pour la prédiction
	 * @param outToclassify: nom du fichier de sortie à prédire
	 */
	void   writeCrossVal(FreqKmer *freqPredict, int percent, string outToclassify);

	/**
	 * Etablit plusieurs validations croisées entre l'objet courant
	 * et celui passé en paramètre
	 * @param freqPredict: objet pour la prédiction
	 * @param percent: pourcentage pour la prédiction
	 * @param i: nombre de validation croisée à produire
	 * @param id: id pour suffixer les fichiers de sortie
	 * 			  pour la prédiction.
	 */
	void   writeNCrossVal(FreqKmer *freqPredict, int percent, int i,string id);

	/**
	 * Genere un ensemble de fichier pour des
	 * multiples validations croisée, afin d'obtenir une courbe
	 * (non moyennée) de vrais positifs pour la suite.
	 * @param sizeSample: taille de l'échantillonage
	 * @param percent: pourcentage de la prédiction
	 * @param start_win_predict: par quelle taille de fenetre commence la prediction
	 * 							 taille minimale du read
	 * @param end: par quelle taille de fenetre termine la prediction
	 * 			   taille maximale du read
	 * @param step: De conmbien incrémenter le prochain read
	 * @param nCross: nombre de validation croisée à produire pour un read.
	 */
	void   generateWekaData(int sizeSample,int percent,int start_win_predict,int end,int step, int nCross);

	/**
	 * Obtenir le nombre de sequence
	 * pour un jeu de données
	 * @param data_i: indice du jeu de données
	 * @return: le nombre de sequences dans data[data_i]
	 */
	int 	getNbSeqData(int data_i);

	/**
	 * Permet d'avoir un format weka
	 * de l'objet courant
	 * @param out: nom du fichier de sortie sans l'extension
	 */
	void 	writeWeka(string out);

};

#endif /* FREQKMER_H_ */
