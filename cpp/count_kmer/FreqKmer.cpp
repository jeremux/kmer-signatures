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
#include <iostream>
#include <algorithm>
using namespace std;

#define LOADALL 1
/******************************************************************************
 *
 *Constructeur
 ****************************************************************************/


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
FreqKmer::FreqKmer(int win_size,int s,bool list, string file,string patternFile, bool b,string key)
{
	printSwitch(getSwitch(0));
	key_fasta=key;
	indexTaxaInFasta=NULL;
	freqFilled=false;
	pathRoot="null";
	nbSeq=NULL;
	data=NULL;
	freq=NULL;
	mask=NULL;
	patterns=NULL;
	kmerSpace=NULL;
	shift=s;
	nCol=0;
	nSeq=0;
	nLine=0;
	nPattern=0;
	nbFastaFile=0;
	winSize=win_size;
	indexLineDataSeq=NULL;
	indexLineData=NULL;
	listData="null";
	noData=b;
	taxaDataSeq=NULL;
	isList=list;
	pathFastaFile=file;
	pathPattern=patternFile;
	initFromRoot=false;
	initWithJump=true;
	/* initialisation des patterns de kmer */
	initPatterns(patternFile);

	/* si c'est une liste en parametre alors on init à partir de initDataFromListFastaPath */
	if (list)
	{
		initDataFromListFastaPath(pathFastaFile);
		initPathFasta(pathFastaFile);
		listData = pathFastaFile;
	}
	else
	{
		initFromFasta(pathFastaFile);
		pathFasta = new string[1];
		pathFasta[0] = file;
	}
}


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
FreqKmer::FreqKmer(int win_size,bool list, string file,string patternFile, bool b,string key)
{
	printSwitch(getSwitch(0));
	key_fasta=key;
	freqFilled=false;
	mask=NULL;
	nbSeq=NULL;
	data=NULL;
	freq=NULL;
	patterns=NULL;
	kmerSpace=NULL;
	nCol=0;
	nSeq=0;
	nLine=0;
	nPattern=0;
	nbFastaFile=0;
	winSize=win_size;
	noData=b;
	taxaDataSeq=NULL;
	indexTaxaInFasta=NULL;
	isList=list;
	pathFastaFile=file;
	pathPattern=patternFile;
	initFromRoot=false;
	initWithJump=false;
	pathRoot="null";
	if(win_size>0)
	{
		shift=(win_size*20)/100;
		if (shift==0)
			shift=1;
	}
	else
	{
		shift=0;
	}
	indexLineDataSeq=NULL;
	indexLineData=NULL;
	initPatterns(patternFile);
	if (list)
	{
		initDataFromListFastaPath(pathFastaFile);
		initPathFasta(pathFastaFile);
		listData = pathFastaFile;
	}
	else
	{
		initFromFasta(pathFastaFile);
		pathFasta = new string[1];
		pathFasta[0] = file;
		listData = "null";
	}
}


/**
 * Contruscteur le decalage = 20% de win_size
 * @param win_size:   taille de la fenetre
 * @param patternFile chemin du fichier contenant les patterns de kmer
 * @param b:          noData, si true alors les données ne sont chargees que lors de leur utilisation
 * @param pathR: 	  chemin du taxon où établir l'apprentissage
 * @param key:        mot cles pour le comptage de frequence (cox1,cox2,...,genomes)
 */
FreqKmer::FreqKmer(int win_size,string patternFile, bool b,string pathR,string key)
{
	printSwitch(getSwitch(0));
	key_fasta=key;
	pathRoot=pathR;
	freqFilled=false;
	mask=NULL;
	data=NULL;
	freq=NULL;
	nbSeq=NULL;
	patterns=NULL;
	kmerSpace=NULL;
	nCol=0;
	nSeq=0;
	nLine=0;
	nPattern=0;
	nbFastaFile=0;
	winSize=win_size;
	noData=b;
	taxaDataSeq=NULL;
	indexTaxaInFasta=NULL;
	isList=true;
	pathPattern=patternFile;
	initFromRoot=true;
	initWithJump=false;
	if(win_size>0)
	{
		shift=(win_size*20)/100;
		if (shift==0)
			shift=1;
	}
	else
	{
		shift=0;
	}
	getDirTaxonFromPath(pathRoot,nbChildTaxa,pathChildTaxa,idTaxa);
	nbChildTaxa = pathChildTaxa.size();
	indexTaxaInFasta = new int[nbChildTaxa];
	initTabIndexTaxaInFasta(key_fasta);
	writeListFasta();
	pathFastaFile=pathR+"/list_fasta.txt";
	indexLineDataSeq=NULL;
	indexLineData=NULL;
	initPatterns(patternFile);
	initDataFromListFastaPath(pathFastaFile);
	initPathFasta(pathFastaFile);
	listData = pathFastaFile;
}


/**
 * Contruscteur
 * @param win_size:   taille de la fenetre
 * @param s:		  taille du decalage
 * @param patternFile chemin du fichier contenant les patterns de kmer
 * @param b:          noData, si true alors les données ne sont chargees que lors de leur utilisation
 * @param pathR: 	  chemin du taxon où établir l'apprentissage
 * @param key:        mot cles pour le comptage de frequence (cox1,cox2,...,genomes)
 */
FreqKmer::FreqKmer(int win_size,int s,string patternFile, bool b,string pathR,string key)
{
	printSwitch(getSwitch(0));
	key_fasta=key;
	pathRoot=pathR;
	freqFilled=false;
	mask=NULL;
	data=NULL;
	freq=NULL;
	nbSeq=NULL;
	patterns=NULL;
	kmerSpace=NULL;
	nCol=0;
	nSeq=0;
	nLine=0;
	nPattern=0;
	nbFastaFile=0;
	winSize=win_size;
	noData=b;
	taxaDataSeq=NULL;
	indexTaxaInFasta=NULL;
	isList=true;
	pathPattern=patternFile;
	initFromRoot=true;
	initWithJump=false;
	shift=s;
	getDirTaxonFromPath(pathRoot,nbChildTaxa,pathChildTaxa,idTaxa);
	nbChildTaxa = pathChildTaxa.size();
	indexTaxaInFasta = new int[nbChildTaxa];
	initTabIndexTaxaInFasta(key_fasta);
	writeListFasta();
	pathFastaFile=pathR+"/list_fasta.txt";
	indexLineDataSeq=NULL;
	indexLineData=NULL;
	initPatterns(patternFile);
	initDataFromListFastaPath(pathFastaFile);
	initPathFasta(pathFastaFile);
	listData = pathFastaFile;
}
/******************************************************************************
 *
 *Destructeur
 ****************************************************************************/
FreqKmer::~FreqKmer()
{
	if(dataVerbose){
		cerr << "Debut FreqKmer::~FreqKmer()\n ";
		cerr.flush();
	}
	for (unsigned int var = 0; var < (unsigned int) nPattern ; var++)
	{
		delete patterns[var];
	}
	delete[] patterns;
	if(!noData)
	{
		for (unsigned int var=0; var < (unsigned int) nbFastaFile ; var++)
		{
			if(data[var]==NULL)
			{
				cerr << "WARNING in ~FreqKmer, NULL pointer on Data object in noData mode\n";
				cerr.flush();
			} else
			{
				delete data[var];
			}
		}
	}
	if (pathRoot!="")
	{
		delete[] indexTaxaInFasta;
	}
	delete[] data;
	for(int i=0;i<nbFastaFile;i++)
	{
		if(indexLineDataSeq[i]!=NULL && mask[i]!=NULL)
		{
			delete[] indexLineDataSeq[i];
			delete[] mask[i];
		}
	}
	delete[] mask;
	delete[] indexLineDataSeq;
	delete[] indexLineData;
	if(pathFasta!=NULL)
		delete[] pathFasta;
	if(dataVerbose){
		cerr << "Debut delete freq\n ";
		cerr.flush();
	}
	if (freq!=NULL)
	{
		for (int var=0; var < nLine ; var++)
		{
			if(freq[var]!=NULL)
			{
				//				cout << "DELETE freq[" << var << "]\n";
				delete[] freq[var];
				freq[var]=NULL;
			}
		}
		delete[] freq;
	}
	if(dataVerbose){
		cerr << "fin delete freq\n ";
		cerr.flush();
	}
	delete[] kmerSpace;
	delete[] nbSeq;
	if(dataVerbose){
		cerr << "Debut FreqKmer::~FreqKmer()\n ";
		cerr.flush();
	}
}


/**
 * Permet d'initialiser le tableau
 * patterns[] a partir du fichier en parametre
 * @param fichier	chemin du fichier de patterns de kmer
 */
void FreqKmer::initPatterns(string fichier)
{
	int tailleLigne=0;
	string ligne;
	ifstream file(fichier.c_str());
	getline(file,ligne);
	tailleLigne = ligne.length();
	int var=0;
	/* On compte combien de pattern on a */
	while (file)
	{
		/* on s'assure que ce n'est pas une ligne vide */
		if(tailleLigne!=0)
		{
			nPattern++;
		}
		/* lecture de la ligne suivante */
		getline(file,ligne);
		tailleLigne = ligne.length();
	}
	/* seconde lecture du fichier */
	ifstream file2(fichier.c_str());
	/* alloue un tableau pour les patterns,
	 * le nombe de pattern est à présent connu
	 */
	patterns = new Pattern*[nPattern];
	getline(file2,ligne);
	tailleLigne = ligne.length();
	/* On rempli le tableau patterns */
	while (file2)
	{
		/* on s'assure que ce n'est pas une ligne vide */
		if(tailleLigne!=0)
		{
			patterns[var++] = new Pattern(ligne);
		}
		/* lecture de la ligne suivante */
		getline(file2,ligne);
		tailleLigne = ligne.length();
	}
	/* Ayant le nombre de kmer
	 * on peut initialiser la map pour l'espace
	 * des kmers (indice de fin de la colonne)
	 */
	kmerSpace = new int[nPattern];
	/* Pour chaque pattern */
	for(int i=0 ; i<nPattern ; i++)
	{
		/* si c'est le premier */
		if (i==0)
		{
			/* alors l'indice de fin pour le premier kmer = (4^k)-1 */
			kmerSpace[i]=patterns[i]->getAllCombi()-1;
		}
		else
		{
			/* sinon c'est l'indice de fin du kmer précent + le nombre de combi du kmer i courant */
			kmerSpace[i]=kmerSpace[i-1]+patterns[i]->getAllCombi();
		}
		/* le nombre total de colonne */
		nCol += patterns[i]->getAllCombi();
	}
}


/**
 * Permet de se deplacer horizontalement
 * @param i:	index du pattern
 * return:	    l'indice du début de la colonne du i-ème pattern
 */
int FreqKmer::obtainStartColKmer(int i)
{
	/* si c'est le premier kmer */
	if (i==0)
	{
		/* alors l'indice de colonne du kmer d'indice 0 commence a 0 */
		return 0;
	}
	else
	{
		/* sinon il commence directement apres l'indice de fin du kmer i-1 */
		return kmerSpace[i-1]+1;
	}
}


/**
 * Permet de se deplacer horizontalement
 * @param i:	index du pattern
 * return:	    l'indice de fin de la colonne du i-ème pattern
 */
int FreqKmer::obtainEndColKmer(int i)
{
	return kmerSpace[i];
}


/**
 * Permet d'initailiser le jeux de donnes
 * @param fichier:	chemin d'un fichier contenant une lsite de chemin fasta
 * Cette methode permet d'etablir la carte du tableau de frequence, en particulier
 * le nombre de ligne
 */
void FreqKmer::initDataFromListFastaPath(string fichier)
{
	if(dataVerbose)
	{
		cerr << "Debut FreqKmer::initDataFromListFastaPath("<< fichier << ")\n";
		cerr.flush();
	}
	idTaxaFromData.erase(idTaxaFromData.begin(),idTaxaFromData.end());
	string ligne;
	int cpt=0;
	int cpt2=-1;
	int tailleSeq = 0;
	ifstream file(fichier.c_str());
	getline(file,ligne);
	int nbFichier = 0;
	int tailleLigne = ligne.length();
	string taxid;
	/* On lit la liste */
	while (file)
	{
		/* Si ce n'est pas une ligne blanche, vide */
		if(tailleLigne!=0)
		{
			nbFichier++;
		}
		/* on lit la ligne suivante */
		getline(file,ligne);
		tailleLigne = ligne.length();
	}
	/* on a notre nombre de fichier, calculé précédemment */
	nbFastaFile=nbFichier;
	/* On peut alors initialiser la premiere dimension de Data et indexLineData */
	data=new Data*[nbFichier];
	indexLineData = new int[nbFichier];
	/* seconde lecture on va initialiser les donnes */
	ifstream file2(fichier.c_str());
	getline(file2,ligne);
	tailleLigne = ligne.length();
	/* J'init à partir de chaque ligne du fichier */
	indexLineDataSeq = new int*[nbFastaFile];
	nbSeq = new int[nbFastaFile];
	mask  = new bool*[nbFastaFile];
	while (file2)
	{
		//		cout << "Traitement de " << ligne << "\n";
		if(tailleLigne!=0)
		{
			/* au premier tour cpt=0, ensuite 1...*/
			//			cout << "TOTOTOTO\n";
			data[cpt] = new Data();
			//			cerr << "lecture de " << ligne << "\n";
			/* J'initilise ma donnée */
			string tmpligne=ligne;
			taxid = getTaxidFromString(ligne);
			idTaxaFromData.push_back(taxid);
			data[cpt]->initFrom(tmpligne,Fasta);
			indexLineDataSeq[cpt] = new int[data[cpt]->getNtaxa()];
			mask[cpt] = new bool[data[cpt]->getNtaxa()];
			/* On incrémente le nombre de sequence total pour le comptage */
			nSeq += data[cpt]->getNtaxa();
			nbSeq[cpt] = data[cpt]->getNtaxa();
			/* traitement pour determiner le nombre de ligne nLine de la table freq */
			if (winSize>0)
			{
				/* Pour chaque sequence du fichier fasta courant */
				for (int var = 0; var < data[cpt]->getNtaxa(); var++)
				{
					/* Je recupere sa taille */
					tailleSeq = data[cpt]->getPrimarySequence(var).length();
					/* Si la fenetre est plus grande que la sequence, alors on a un seul comptage à faire */
					if (tailleSeq<winSize)
					{
						nLine += 1;
					}
					else
					{
						nLine += obtainNbLineWindow(0,tailleSeq-1,winSize,shift);
					}
					/* Je peux à cette étape savoir l'indice de fin de ligne pour
					 * la var-ème sequence du cpt-ème jeu de données
					 */
					indexLineDataSeq[cpt][var]=nLine-1;
					mask[cpt][var]=true;
				}
			}
			/* Si winSize <= 0 alors la taille de la fenetre est celle de la sequence */
			else
			{
				nLine +=  data[cpt]->getNtaxa();
				for(int i=0;i<data[cpt]->getNtaxa();i++)
				{
					indexLineDataSeq[cpt][i] = cpt2++;
					mask[cpt][i]=true;
				}
			}
			indexLineData[cpt]=nLine-1;
			if(noData)
			{
				// cerr << "Fermeture Data\n";
				delete data[cpt];
				data[cpt]=NULL;
				//				cout << "data[" << cpt << "] = NULL \n";
			}
			cpt++;
		}
		getline(file2,ligne);
		tailleLigne = ligne.length();
	}
	if(dataVerbose)
	{
		cerr << "Fin FreqKmer::initDataFromListFastaPath("<< fichier << ")\n";
		cerr.flush();
	}
}


/** Permet de garder en memoire la liste des chemins fasata
 * @param fichier:	chemin du fichier de la liste contenant les chemins fasta
 * Cette methode permet de recuperer les taxids de chaque chemins
 */
void FreqKmer::initPathFasta(string fichier)
{
	string ligne;
	string taxid;
	int cpt=0;
	ifstream file(fichier.c_str());
	idTaxaFromData.erase(idTaxaFromData.begin(),idTaxaFromData.end());
	getline(file,ligne);
	int tailleLigne = ligne.length();
	/* On lit la liste */
	pathFasta = new string[nbFastaFile];
	while (file)
	{
		/* Si ce n'est pas une ligne blanche, vide */
		if(tailleLigne!=0)
		{
			taxid = getTaxidFromString(ligne);
			idTaxaFromData.push_back(taxid);
			pathFasta[cpt]=ligne;
		}
		/* on lit la ligne suivante */
		getline(file,ligne);
		tailleLigne = ligne.length();
		cpt++;
	}
}


/**
 * Permet d'initailiser le jeux de donnes
 * @param fichier:	fichier fasta
 * Cette methode permet d'etablir la carte du tableau de frequence, en particulier
 * le nombre de ligne
 */
void FreqKmer::initFromFasta(string fichier)
{
	nbFastaFile=1;
	data=new Data*[1];
	indexLineData = new int[1];
	int cpt2=-1;
	data[0] = new Data();
	data[0]->initFrom(fichier,Fasta);
	nSeq = data[0]->getNtaxa();
	indexLineDataSeq = new int *[1];
	indexLineDataSeq[0] = new int[nSeq];
	mask  = new bool *[1];
	nbSeq = new int[1];
	nbSeq[0] = nSeq;
	mask[0] = new bool[nSeq];
	int taille=0;
	/* si une taille de fenetre est définie */
	if (winSize>0)
	{
		/* alors on peut calculer le nombre de ligne */
		/* Pour chaque sesquence du fasta courant */
		for (int var = 0; var < nSeq; var++)
		{
			taille= data[0]->getPrimarySequence(var).length();
			if (taille<winSize)
			{
				nLine += 1;
			}
			else
			{
				//nLine += taille-winSize+1;
				nLine += obtainNbLineWindow(0,taille-1,winSize,shift);
			}
			indexLineDataSeq[0][var]=nLine-1;
			mask[0][var]=true;
		}
	}
	else
	{
		for (int var = 0; var < nSeq; var++)
		{
			indexLineDataSeq[0][var]=cpt2++;
			mask[0][var]=true;
		}
		nLine = nSeq;
	}
	indexLineData[0] = nLine-1;
	if(noData)
	{
		delete data[0];
		data[0]=NULL;
	}
}


/**
 * Methode permettant de savoir à quelle colonne
 * correspond un kmer donnée dans une sequence
 * @param indexPattern:		indice du pattern dans le tableau patterns
 * @param seq:				sequence pour laquelle on effectue le comptage
 * @param pos:				indice où l'on se trouve dans la sequence
 * return: 					indice de la colonne dans la table freq
 */
int FreqKmer::obtainColIndex(int indexPattern,int *seq,int pos)
{
	return obtainStartColKmer(indexPattern)+patterns[indexPattern]->getKmer(seq,pos);
}


/**
 * initialise la première dimension du tableau freq.
 */
void FreqKmer::initFreq()
{
	if(dataVerbose){
		cerr << "Debut FreqKmer::initFreq()\n ";
		cerr.flush();
	}
	freq = new double *[nLine];
	for(int i=0 ; i<nLine ; i++)
	{
		freq[i] = NULL;
	}
	if(dataVerbose){
		cerr << "Fin FreqKmer::fillFreq()\n ";
		cerr.flush();
	}
}


/***
 * effectue le comptage dans une sous sequence
 * @param	seq: sequence où compter
 * @param	win_length: taille de la sous sequence
 * @param	pos: où commencer à compter
 * @param	indexPattern: indice du kmer courant
 * @param	previous: buffer contenant la frequence en kmer de
 * 					  la partie commune entre la sequence courante et precedent
 * @param   start_line: à partir de quel ligne ecrire.
 */
void FreqKmer::winCount(int *seq,int win_length,int pos ,int indexPattern,int *previous,int start_line)
{
	if(dataVerbose){
		cerr << "Debut FreqKmer::winCount\n ";
		cerr << " win_length = " << win_length << "\t pos = " << pos << "\n";
		cerr.flush();
	}
	int col;
	/* on s'arrête à la dernière fenetre possible */
	int fin = pos+win_length-patterns[indexPattern]->getSizePattern();
	for(int i = pos; i <= fin ; i++)
	{
		/* indice du kmer à la position i */
		col = obtainColIndex(indexPattern,seq,i);
		if(col >= 0 && col < getNCol())
		{
			/* si ligne non allouee */
			if(freq[start_line]==NULL)
			{
				/* allocation */
				freq[start_line] = new double[nCol];
				for(int p=0 ; p<nCol ; p++)
				{
					freq[start_line][p]=0.0;
				}
			}
			freq[start_line][col]+=1;
		}
		/* Je sauvegarde le premier comptage dans un buffer */
		/* Au premier tour i-pos=0, ensuite 1... */
		previous[i-pos]=col;
	}
	if(dataVerbose)
	{
		cerr << "Fin FreqKmer::winCount\n ";
		cerr.flush();
	}
}


/**
 * Permet de copier un la partie commune
 * d'un buffer vers un autre et de mettre a jour le nouveau buffer
 * @param previous:		buffer du comptage de la fenetre precedente
 * @param current:		buffer a mettre a jour
 * @param buf_size:		taille des buffer
 * @param indexPattern:	index du pattern courant dans patterns
 * @param seq: 			sequence pour laquelle on effectue le comptage
 * @param pos:			position ou l'on se trouve dans la sequence
 */
void FreqKmer::copyBuffAndCount(int *current,int *previous,int buf_size, int indexPattern,int *seq,int pos)
{
	for(int i=0;i<buf_size;i++)
	{
		current[i]=-1;
	}
	int cpt=0;
	/* On commence au décalage et on s'arrête à la fin
	 * du buffer previous
	 */
	for(int i=shift;i<buf_size;i++)
	{
		/* Au depart i-shift = 0, ensuite 1...*/
		current[i-shift]=previous[i];
		cpt++;
	}
	/* On reprend là où
	 * on s'est arrêté dans la boucle précédente
	 * c'est à dire i=buf_size, et on rempli les
	 * shift case restante.
	 */
	// si le buffer est plus petit que le
	// 	decalage on doit tout recompter
	if (buf_size < shift)
	{
		for(int i=0;i<buf_size;i++)
		{
			current[i] = obtainColIndex(indexPattern,seq,pos+i);
		}
	}
	else
	{
		for(int i=pos+buf_size-shift;i<=pos+winSize-patterns[indexPattern]->getSizePattern();i++)
		{
			current[cpt++] = obtainColIndex(indexPattern,seq,i);
		}
	}
}


/**
 * Methode permettant d'echanger deux buffer
 * @param current, previous: les buffers a intervertir
 * @param buf_size:			 taille des buffers
 */
void FreqKmer::swap(int *current,int *previous,int buf_size)
{
	int *tmp = new int[buf_size];
	for(int i=0;i<buf_size;i++)
	{
		tmp[i]=previous[i];
		previous[i]=current[i];
		current[i]=tmp[i];
	}
	delete[] tmp;
	tmp=NULL;
}


/**
 * Methode implementee
 * pour le debogage, elle permet
 * d'afficher sur la sortie derreur le contenu
 * d'un buffer
 * @param buf:		buffer a afficher
 * @param buf_size:	taille du buffer
 */
void FreqKmer::printBuf(int *buf,int buf_size)
{
	cerr << "[indice - val]\n";
	cerr << "[";
	for(int i=0;i<buf_size-1;i++)
	{
		cerr << i <<" - "<< buf[i] << "][";
	}
	cerr << buf_size-1 <<" - "<<buf[buf_size-1] << "]\n";
}


/**
 * Effectue le comptaga pour seq
 * @param 	seq: 			la sequence ou compter
 * @param 	seq_length : 	taille de la sequence
 * @parm	indexPattern : 	indice du kmer actuel
 * @param 	start_line:		indice de la ligne dans le tableau freq
 */
void FreqKmer::count(int *seq,int seq_length,int indexPattern,int start_line)
{
	if(dataVerbose)
	{
		cerr << "FreqKmer::count()\n ";
		cerr.flush();
	}
	int i=0;
	int j=0;
	int *previous = NULL;
	int *current = NULL;
	int index = start_line;
	/*taille de la sous sequence où travailler */
	int win_length = 0;
	int z;
	/* Initialisation
	 * des tableaux premier
	 * et dernier
	 */
	if (winSize==-1 || winSize>seq_length)
	{
		win_length = seq_length;
	}
	else
	{
		win_length = winSize;
	}
	if(dataVerbose)
	{
		cerr << "--------alloc buffer\n ";
		cerr.flush();
	}
	int buf_size = win_length-patterns[indexPattern]->getSizePattern()+1;
	previous = new int[buf_size];
	current = new int[buf_size];
	/* on recupère le nombre de ligne pour la sequence
	 * afin d'iterer le bon nombre de fois
	 */
	z=obtainNbLineWindow(0,seq_length-1,win_length,shift);
	while(j < z)
	{
		/* on compte seulement si la ligne n'a pas été calculée */
		/* Si j > 0 alors on a déjà effectué le comtpage pour la première fenetre
		 * on utilise l'astuce de decalage
		 */
		if (j>0)
		{
			/* On copie la bonne partie du buffer et on compte les nouveaux kmers */
			copyBuffAndCount(current,previous,buf_size,indexPattern,seq,i);
			/* on increment les kmers trouvés */
			for(int i=0;i<buf_size;i++)
			{
				/* On ne compte que les acgt */
				if(current[i]<getNCol() && current[i]>0)
				{
					/* on instancie la ligne */
					if(freq[index]==NULL)
					{
						//							cout << "ALLOC freq[" << index << "]\n";
						freq[index] = new double[nCol];
						for(int p=0 ; p<nCol ; p++)
						{
							freq[index][p]=0.0;
						}
					}
					freq[index][current[i]]+=1;
				}
			}
			/* le courant devient le precedent */
			swap(current,previous,buf_size);
		}
		else
		{
			//			cout << "first window\n";
			//			cout.flush();
			winCount(seq,win_length,i,indexPattern,previous,start_line);
			//			cout << "end first window\n";
			//			cout.flush();
		}
		/* On decale de shift nucleotides */
		i = i + shift;
		j = j + 1;
		index++;
	}
	if (previous)
	{
		delete[] previous;
		previous = NULL;
	}
	if (current)
	{
		delete[] current;
		current = NULL;
	}
}


/**
 * Rempli le tableau de frequence
 */
void FreqKmer::fillFreq()
{
	if(dataVerbose){
		cerr << "FreqKmer::fillFreq()\n ";
		cerr.flush();
	}
	int nTaxa;
	if(freq==NULL)
		initFreq();
	/* Pour chaque fichier fasta */
	for(int i=0;i<nbFastaFile;i++)
	{

		nTaxa = nbSeq[i];

		for(int j=0;j<nTaxa;j++)
		{
			/* remplir pour la j-ème sequence du i-ème data */
			fillFreq(i,j);
		}
	}
	freqFilled=true;
}


/**
 * Permet d'avoir le nombre de ligne (decalage)
 * pour une sous sequence selon une taille de fenetre
 * @param	i	indice du debut (varie de 0 à taille de la sequence  -1 )
 * @param	j 	indice de fin
 * @param	l	taille de la fenetre
 * @return	le nombre de ligne.
 */
int FreqKmer::obtainNbLineWindow(int i,int j,int l,int pas)
{
	int res = 0;
	if(pas!=0)
	{
		res = ((j-i+1-l)/pas)+1;
	}
	else
	{
		res = 1;
	}
	return res;
}


/**
 * Permet d'avoir le nombre de ligne
 * pour la donnée Data à l'indice i
 * (soit un ensemble de sequence)
 * @param i	indice du data à considérer
 * @return le nombre de ligne pour data[i]
 */
int FreqKmer::obtainNbLineData(int i)
{
	if(i==0)
	{
		return indexLineData[i]+1;
	}
	else
	{
		return indexLineData[i]-indexLineData[i-1];
	}
}


/**
 * Permet d'avoir le nombre de ligne
 * pour la donnée Data à l'indice i
 * (soit un ensemble de sequence)
 * @param i	indice du data à considérer
 * @return le nombre de ligne pour data[i]
 */
int FreqKmer::obtainStartLineData(int i)
{
	if(i==0)
	{
		return 0;
	}
	else
	{
		return indexLineData[i-1]+1;
	}
}

/**
 * Methode pour obtenir l'indice de fin
 * de la ligne pour un data donnée
 * @param i:	i-ème data
 * return:		l'indice de ligne de fin dans la table freq de data[i]
 */
int FreqKmer::obtainEndLineData(int i)
{
	return indexLineData[i];
}

/**
 * Methode pour obtenir l'indice de
 * debut pour une séquence d'un jeu de donnée
 * @param i:	i-ème jeu de donnee
 * @param j:	j-ème sequence
 * @return:		l'indice de ligne de debut dans la table freq de la j-ème sequence de data[i]
 */
int FreqKmer::obtainStartLineDataSeq(int i,int j)
{
	int tmp = 0;
	int tmp2 = i-1;
	if(j==0)
	{
		if(i==0)
		{
			return 0;
		}
		else
		{

			tmp = nbSeq[tmp2]-1;

			return indexLineDataSeq[i-1][tmp]+1;
		}
	}
	else
	{
		return indexLineDataSeq[i][j-1]+1;
	}
	return 0;
}

/**
 * Methode pour obtenir l'indice de
 * fin pour une séquence d'un jeu de donnée
 * @param i:	i-ème jeu de donnee
 * @param j:	j-ème sequence
 * @return:		l'indice de ligne de fin dans la table freq de la j-ème sequence de data[i]
 */
int FreqKmer::obtainEndLineDataSeq(int i,int j)
{
	return indexLineDataSeq[i][j];
}

/**
 * Methode pour obtenir lenombre de ligne
 * pour une séquence d'un jeu de donnée
 * @param i:	i-ème jeu de donnee
 * @param j:	j-ème sequence
 * @return:		le nombre de ligne q'occupe la j-ème sequence de data[i] dans la table freq
 */
int FreqKmer::obtainNbLineDataSeq(int i,int j)
{
	return obtainEndLineDataSeq(i,j)-obtainStartLineDataSeq(i,j)+1;
}

/**
 * Methode pour recuperer l'indice de la sequence et du jeu de
 * donnee à partir du ligne de la table freq
 * @param line:		indice de ligne dans la table freq
 * @param idData:	où sauvegarder l'indice du jeu de donnee
 * @param idSeq:	où sauvegarder l'indice de la sequence
 */
void FreqKmer::obtainDataSeqFromLine(int line,int &idData,int &idSeq)
{
	bool res = false;
	int n;
	for(int i=0;i<nbFastaFile;i++)
	{
		n = nbSeq[i];
		for(int j=0;j<n;j++)
		{
			if(line >= obtainStartLineDataSeq(i,j) && line <= obtainEndLineDataSeq(i,j))
			{
				idData = i;
				idSeq = j;
				res = true;
			}
			if (res)
				break;
		}

		if(res)
			break;
	}
}

/**
 * Methode pour recuperer l'indice de la sequence et du jeu de
 * donnee à partir du ligne de la table freq et egalement la fenetre
 * @param line:		indice de ligne dans la table freq
 * @param idData:	où sauvegarder l'indice du jeu de donnee
 * @param idSeq:	où sauvegarder l'indice de la sequence
 * @param idWin: 	où sauvegarder l'indice de la fenetre
 */
void FreqKmer::obtainDataSeqWinFromLine(int line,int &idData,int &idSeq,int &idWin)
{
	bool res = false;
	int n;
	for(int i=0;i<nbFastaFile;i++)
	{
		n = nbSeq[i];
		for(int j=0;j<n;j++)
		{
			if(line >= obtainStartLineDataSeq(i,j) && line <= obtainEndLineDataSeq(i,j))
			{
				idData = i;
				idSeq = j;
				idWin = line-obtainStartLineDataSeq(i,j)+1;
				res = true;
			}
			if (res)
				break;
		}

		if(res)
			break;
	}
}

/**
 * Methode pour savoir si on
 * considere une sequence dans nos calcul
 * @param i:	indice du jeu de donnee a questionner
 * @param j:	indice de la sequence du i-ème jeu de donnes
 */
bool FreqKmer::takeDataSeq(int i,int j)
{
	return mask[i][j];
}

/**
 * Methode permettant de "cacher" une sequence
 * @param i:	indice du jeu de donnee a cacher
 * @param j:	indice de la sequence du i-ème jeu de donnes
 */
void FreqKmer::setFalseMask(int i,int j)
{
	mask[i][j] = false;
}

/**
 * Methode permettant de recuperer la liste
 * des sous taxons du taxon courant
 * @param dir:		chemin du taxon courant
 * @param nbChild:	où sauvegarder le nombre de sous taxon
 * @param files:	liste de string ou sera enregistré les chemins vers les fastas
 * @param taxids:	liste de string contenant les taxids des sous taxons du taxon courant
 */
int FreqKmer::getDirTaxonFromPath(string dir,int &nbChildTaxa, vector<string> &files,vector<string> &taxids)
{
	files.erase(files.begin(),files.end());
	taxids.erase(taxids.begin(),taxids.end());
	DIR *dp;
	string s1 = "";
	int x = -1;
	unsigned found;
	struct dirent *dirp;
	if((dp  = opendir(dir.c_str())) == NULL)
	{
		cout << "Error(" << errno << ") opening " << dir << endl;
		return errno;
	}

	/* Pour chaque dossier dirp du dossier courrant*/
	while ((dirp = readdir(dp)) != NULL)
	{
		s1 = string(dirp->d_name);
		/* si le nom de ce dossier contient "__" ou est others */
		if (s1.find("__") != std::string::npos || s1=="others")
		{
			/* on push le dossier (après concatenation avec le chemin du dossier courant */
			files.push_back(dir+"/"+s1);
		}
	}

	/* on trie par ordre alphabetique la liste des fichiers fasta
	 * puisque l'ordre n'est pas garantie par l'appel readdir, cela
	 * permet d'avoir des tests valides sur toutes machines.
	 * A noter que le trie n'est pas nécéssaire au fonctionnement du programme.
	 */
	sort(files.begin(),files.end());

	/* Après avoir trié, on recupere les taxids */
	for(unsigned i=0;i< files.size();i++)
	{
		s1 = files[i] ;
		/* on cherche l'indice du dernier slash dans le path */
		found=s1.find_last_of("/\\");

		/* si c'est un slash en fin de path */
		if(found==s1.size()-1)
		{
			// on enlève le slash si en fin de path : ex taxon__A/taxon__B/taxon__C/ "
			// ensuite vaut taxon__A/taxon__B/taxon__C
			s1 = s1.substr(found,1);
			found=s1.find_last_of("/\\");
		}

		/* on recupere le nom du taxon dans l'exemple
		 * ça serait taxon__C
		 */
		s1=s1.substr(found+1,s1.size());
		/* si c'est le dossier others alors pas de taxid */
		if (s1=="others")
		{
			taxids.push_back(s1);
		}
		/* sinon on recupere le taxid */
		else
		{
			/* dans l'exemple de s1 = taxon__C
			 * s1.find("__") = 5
			 * x = 7
			 */
			x = s1.find( "__" ) + 2 ;

			/* Dans l'exemple s1.substr(x,s1.length()) = C */
			taxids.push_back(s1.substr(x,s1.length()));
		}
	}

	/* le nombre de dossier taxon */
	nbChildTaxa = files.size();
	closedir(dp);
	return 0;
}

/**
 * Permet de recuperer le taxid
 * à partir du dossier d'indice dir_i
 * @param dir_i:	indice du dossier
 * return:			le taxid du dossier d'indice dir_i
 */
string FreqKmer::getIdTaxa(int dir_i)
{
	return idTaxa[dir_i];
}

/**
 * Permet de recupere le chemin du
 * dossier d'indice dir_i
 * @param dir_i:	indice du dossier
 * return: 			le chemin du dossier d'indice dir_i
 */
string FreqKmer::getPathChildTaxa(int dir_i)
{
	return pathChildTaxa[dir_i];
}

/**
 * Permet d'avoir le numéro de ligne (qui
 * commence à 0) de début
 * dans la liste des chemin fasta pour un taxon donnée
 * @param i:	index du taxon
 * return:		l'indice de début de ligne du taxon i
 */
int FreqKmer::obtainStartLineTaxaInFastaList(int i)
{
	if(i==0)
	{
		return 0;
	}
	else
	{
		return indexTaxaInFasta[i-1]+1;
	}
}

/**
 * Permet d'avoir le numéro de ligne (qui
 * commence à 0) de fin
 * dans la liste des chemin fasta pour un taxon donnée
 * @param i:	index du taxon
 * return:		l'indice de fin de ligne du taxon i
 */
int FreqKmer::obtainEndLineTaxaInFastaList(int i)
{
	return indexTaxaInFasta[i];
}

/**
 * Permet d'avoir le nombre de fichier fasta
 * pour un taxon donnée
 * @param i:	index du taxon
 * return:		le nombre de fichier fasta pour le i-ème taxon
 */
int FreqKmer::obtainNbLineTaxaInFastaList(int i)
{
	int res = 0;
	res = (obtainEndLineTaxaInFastaList(i)-obtainStartLineTaxaInFastaList(i))+1;
	return res;
}

/**
 * Initialise la liste de fichier fasta
 * des sous taxons du taxon (dossier) courant
 * @param key_fasta:	mot-clé de la sequence a recuperer
 * 						cox1,cox2,...,genomes.
 */
void FreqKmer::initTabIndexTaxaInFasta(string key_fasta)
{
	struct dirent *dirp;
	DIR *dp;
	int cpt=-1;
	string s1="";
	string path="";
	string path_tmp = "";
	string file_name = "";
	string extension = "";
	string ancien = "";
	int index=-1;
	string current="";
	string previous="";
	string taxid;
	unsigned found;
	int x = -1;
	//	cerr << "nbChildttaxa = " << nbChildTaxa << "\n";
	listPathFasta.erase(listPathFasta.begin(),listPathFasta.end());

	/* Pour chaque sous taxon du taxon courant */
	for(int i=0; i < nbChildTaxa ; i++)
	{
		/* on recupere le chemin du taxon courant */
		path = getPathChildTaxa(i);

		/* dossier ou se trouve les sequences à recuperer */
		path += "/data/fasta/nucleotides/" + key_fasta;

		if((dp  = opendir(path.c_str()) )== NULL)
		{
			cerr << "Error(" << errno << ") opening " << path << endl;
			exit(0);
		}

		/* Pour chaque fichier du dossier path */
		while ((dirp = readdir(dp)) != NULL)
		{
			file_name=dirp->d_name;
			path_tmp = path+"/"+file_name;

			/* Si ce n'est pas un dossier */
			if(!(dirp->d_type == DT_DIR))
			{
				/* On recupere l'extension */
				x = file_name.find( "." ) + 1 ;
				extension = file_name.substr(x,file_name.length());
				/* si c'est un fasta on push */
				if(extension=="fasta")
				{
					listPathFasta.push_back(path_tmp);
					cpt++;
				}
			}
		}

		closedir(dp);
	}

	/**
	 * On trie la list listPathFasta car l'appel
	 * readdir ne garantie pas d'ordre et on souhaite que
	 * les tests soient indépendantes des machines
	 */
	sort(listPathFasta.begin(),listPathFasta.end());

	/* un marqueur sur l'ancien, au début egale au premier fasta */
	/* ex: ancien = taxon__A/taxon__B/cox1-xxx.fasta */
	ancien=listPathFasta[0];


	found=ancien.find_last_of("/\\");

	/* dans l'ex previous = taxon__A/taxon__B */
	previous=ancien.substr(0,found);

	/**
	 * On fait la meme chose pour les autres
	 * chemins et si on tombe sur un chemin different
	 * on peut mettre a jour la table indexTaxaInFasta
	 */
	for(unsigned int j=1;j<listPathFasta.size();j++)
	{
		current=listPathFasta[j];
		found=current.find_last_of("/\\");
		current=current.substr(0,found);
		if(current!=previous)
		{
			previous=current;
			indexTaxaInFasta[++index]=j-1;
		}
	}
	indexTaxaInFasta[++index]=listPathFasta.size()-1;
}

/**
 * Methode permettant d'ecrire en
 * dur dans le dossier courant la liste des fasta
 * des sous taxon du dossier courant
 */
void FreqKmer::writeListFasta()
{
	if (dataVerbose)
	{
		cerr << "Debut FreqKmer::writeListFasta()\n";
	}
	string output = pathRoot+"/"+ "list_fasta.txt";
	ofstream myfile ;
	myfile.open(output.c_str());
	for(unsigned int j=0;j<listPathFasta.size();j++)
	{

		myfile << listPathFasta[j] << "\n";
	}
	myfile.close();
	if (dataVerbose)
	{
		cerr << "Fin FreqKmer::writeListFasta()\n";
	}
}

/**
 * Tirage aleatoire sans remise
 * @param result:		liste des entiers tirés
 * @param tabSize:		borne sup du tirage
 * @param sampleSize:	nombre d'entier à tirer
 */
void FreqKmer::randomTab(vector<int> *result,int tabSize,int sampleSize)
{
	if(dataVerbose)
	{
		cerr << "DEBUT FreqKmer::randomTab(tab,"<<tabSize<<","<<sampleSize<<")\n";
	}
	result->erase(result->begin(),result->end());

	if(sampleSize>tabSize)
	{
		for(int i=0;i<tabSize;i++)
		{
			result->push_back(i);
		}

		return;
	}

	int *tmp = new int[tabSize];
	int r = -1;

	/* borne sup du tirage */
	int sup = tabSize;
	int val_tmp;

	/*init de chaque case */
	for(int i=0;i<tabSize;i++)
	{
		tmp[i]=i;
	}

	/* on fait sampleSize tirage */
	/* invariant:
	 * dans le tableau tmp
	 * entre 0 et sup-1 on a des entiers
	 * non tirés
	 */
	for(int j=0;j<sampleSize;j++)
	{
		/* au debut on tire entre
		 * 0 et sup=tabSize exclu, puis entre
		 * 0 et sup-1
		 * 0 et sup-2
		 * ...
		 */
		r = rand() % sup;
		val_tmp = tmp[sup-1];
		/* on recupère l'entier tiré */

		result->push_back(tmp[r]);
		/* on range l'entier tiré à la fin du tableau tmp */
		tmp[sup-1] = tmp[r];
		tmp[r] = val_tmp;
		sup--;
	}
	delete[] tmp;
	if(dataVerbose)
	{
		cerr << "FIN FreqKmer::randomTab(tab,"<<tabSize<<","<<sampleSize<<")\n";
	}
}

/**
 * Permet d'obtenier un echantillon du jeu
 * de donnees initialement fourni
 * @param sampleSize:	taille de l'échantillon
 * return:				un pointeur d'objet de type FreqKmer
 */
FreqKmer* FreqKmer::sampleMe(int sampleSize)
{
	FreqKmer *res;
	bool b = noData;
	int n;
	sampledTaxon.erase(sampledTaxon.begin(),sampledTaxon.end());

	if (dataVerbose)
	{
		cerr << "Debut FreqKmer::sampleMe("<< sampleSize << ")\n";
	}

	if(initFromRoot)
	{
		if(initWithJump)
		{
			res = new FreqKmer(winSize,shift,pathPattern,b,pathRoot,key_fasta);
		}
		else
		{
			res = new FreqKmer(winSize,pathPattern,b,pathRoot,key_fasta);
		}
	}
	else
	{
		if(initWithJump)
		{
			res = new FreqKmer(winSize,shift,isList,pathFastaFile,pathPattern,b,key_fasta);
		}
		else
		{
			res = new FreqKmer(winSize,isList,pathFastaFile,pathPattern,b,key_fasta);
		}
	}

	if(LOADALL && noData)
	{
		cerr << "Changement en cours de route noData = false\n";
		res->setNoData(false);
		res->loadAll();
	}

	int nbSequences = 0;
	int nbSeqTaxa = 0;
	int d,f;
	d=0;
	f=0;
	bool **mask_tmp;
	res->freqFilled=this->freqFilled;

	/* on recupére les lignes calculées */
	if(this->freqFilled)
	{
		res->initFreq();
		for(int i=0;i<nLine;i++)
		{
			if(res->freq[i]==NULL)
			{
				res->freq[i] = new double[nCol];
				for(int p=0 ; p<res->nCol ; p++)
				{
					res->freq[i][p]=0.0;
				}
			}
			for(int j=0;j<nCol;j++)
			{
				res->freq[i][j]=this->freq[i][j];
			}
		}
	}
	mask_tmp = new bool*[nbFastaFile];
	vector<int> candidates;
	for(int i=0;i<nbFastaFile;i++)
	{

		nbSequences = nbSeq[i];
		mask_tmp[i] = new bool[nbSequences];

		/* on initialise le mask temporaire à false */
		for(int j=0;j<nbSequences;j++)
		{
			mask_tmp[i][j]=false;

		}


	}

	for(int i=0;i<nbChildTaxa;i++)
	{
		nbSeqTaxa = getNSeqInTaxa(i);
		/* si on doit tirer plus de qu'il y a
		 * de seq alors on tire tout
		 */

		if(sampleSize>=nbSeqTaxa || sampleSize==-1)
		{
			d = obtainStartLineTaxaInFastaList(i);
			f = obtainEndLineTaxaInFastaList(i);
			for(int j=d; j<=f;j++)
			{
				n = nbSeq[j];

				for(int k=0;k<n;k++)
				{
					mask_tmp[j][k]=true;
					/* on garde en memoire les donnees recuperees */
					sampledTaxon.push_back(intPair(j,k));
				}

			}
		}
		/* on peut alors tirer n seq dans le taxon i */
		else
		{
			/* on tire sampleSize parmi les nbSeqTaxa sequences
			 * dans le taxon courant
			 * (tirage sans remise pour les seq.!)
			 */
			randomTab(&candidates,nbSeqTaxa,sampleSize);
			/* je mets à jour mon mask selon les candidats trouvés */
			maskTab(&candidates,mask_tmp,i);
		}
	}
	/* MAJ de la table mask */
	for(int i=0;i<nbFastaFile;i++)
	{

		nbSequences = nbSeq[i];

		for(int j=0;j<nbSequences;j++)
		{
			res->mask[i][j]=mask_tmp[i][j];
			/* on garde une trace des sequences gardées, donc quand mask_tmp[i][j]=true */
			if(mask_tmp[i][j])
			{
				sampledTaxon.push_back(intPair(i,j));
			}
		}
		delete[] mask_tmp[i];
	}
	delete[] mask_tmp;
	if (dataVerbose)
	{
		cerr << "Fin FreqKmer::sampleMe("<< sampleSize << ")\n";
	}
	/* on retourne notre nouvel objet */
	return res;
}


/**
 * Permet d'avoir le nombre de sequence
 * d'un taxon
 * @param i:	indice du taxon
 * return:		le nombre de sequence du i-ème sous taxon du taxon courant
 */
int FreqKmer::getNSeqInTaxa(int i)
{
	int nbData = 0;

	/* nombre de fichier fasta */
	nbData = obtainNbLineTaxaInFastaList(i) ;
	//	cerr << "nbData = " << nbData << "\n";
	int res=0;
	int startLineInFasta = 0;

	/* indice de début dans data */
	startLineInFasta = obtainStartLineTaxaInFastaList(i);

	/* Pour chaque fichier fasta */
	for(int j=0; j<nbData ; j++)
	{
		res+=nbSeq[j+startLineInFasta];
	}
	return res;
}

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
int FreqKmer::getNbTrue(bool *tab,int start,int end)
{
	int res=0;
	for(int i=start;i<=end;i++)
	{
		if(tab[i])
		{
			res++;
		}
	}
	return res;
}

/**
 * Methode permettant de savoir combien
 * de sequences sont considérées, n'est
 * pas utile au programme mais permet des tests
 * basiques
 */
int FreqKmer::getNbAllTrue()
{
	int res = 0;
	for(int i=0;i<nbFastaFile;i++)
	{

		res+=getNbTrue(mask[i],0,nbSeq[i]-1);

	}
	return res;
}

/**
 *
 */
void FreqKmer::maskTab(vector<int> *candidate,bool **mask_tmp, int indexTaxa)
{
	/*TODO: no data */
	int d = obtainStartLineTaxaInFastaList(indexTaxa);
	int f = obtainEndLineTaxaInFastaList(indexTaxa);
	int cpt=0;
	int nSeq = 0;
	for(int i=d; i <= f ; i++)
	{
		/* TODO: no data instancier data ici */

		nSeq = nbSeq[i];

		for(int j=0 ; j < nSeq ; j++)
		{
			if(find(candidate->begin(), candidate->end(),cpt) != candidate->end())
			{

				/* data[i][j] is a candidate for sample */
				mask_tmp[i][j]=true;
			}
			cpt++;
		}
	}
}


double FreqKmer::getSum(int indexPattern,int indexLine)
{
	double s=0;
	//cout << "calcul pour i = " << obtainStartColKmer(indexPattern) << " à " << obtainEndColKmer(indexPattern) << " : ";
	if(freq[indexLine]!=NULL)
	{
		for(int i=obtainStartColKmer(indexPattern);i<=obtainEndColKmer(indexPattern);i++)
		{
			s += freq[indexLine][i];
		}
	}
	return s;
}
void FreqKmer::normalizeLine(int indexPattern,int line)
{
	double s = getSum(indexPattern,line);
	int nbComb ;
	if(normalNormalize)
	{
		nbComb=1;
	}
	else
	{
		nbComb = patterns[indexPattern]->getAllCombi();
	}


	if (s!=0)
	{
		for(int i=obtainStartColKmer(indexPattern);i<=obtainEndColKmer(indexPattern);i++)
		{
			if(freq[line]!=NULL)
				freq[line][i] = (nbComb) * freq[line][i]/s;
		}
	}
}
void FreqKmer::normalize()
{
	if(freqFilled)
	{
		for(int i=0;i<nPattern;i++)
		{
			for(int j=0;j<nLine;j++)
			{
				if(freq[j]!=NULL)
					normalizeLine(i,j);
			}
		}
	}
}
void FreqKmer::writeConfFeq(string output)
{
	ofstream myfile ;
	myfile.open(output.c_str());
	//	myfile << "INIT_CONS\n";
	if(initFromRoot)
	{
		if(initWithJump)
		{
			myfile << "4\n" << winSize << "\n"<< shift << "\n";
			myfile << pathPattern << "\n" << noData << "\n" << pathRoot << "\n" << key_fasta << "\n";
		}
		else
		{
			myfile << "3\n" << winSize <<  "\n";
			myfile << pathPattern << "\n" << noData << "\n" << pathRoot << "\n" << key_fasta << "\n";
		}
	}
	else
	{
		if(!initWithJump)
		{
			myfile << "2\n" << winSize << "\n"<< isList << "\n" << pathFastaFile << "\n" ;
			myfile << pathPattern << "\n" << noData << "\n"  << key_fasta << "\n";
		}
		else
		{
			myfile << "1\n" << winSize << "\n" << shift << "\n"<< isList << "\n" << pathFastaFile << "\n" ;
			myfile << pathPattern << "\n" << noData << "\n"  << key_fasta << "\n";
		}
	}
	myfile << "FIN_CONS\n";
	//	myfile << "DIM\n";
	myfile << nLine << "\n";
	myfile << nCol << "\n";
	myfile << "FIN_DIM\n";
	//	cout << endl;
	for(int i=0;i<nLine;i++)
	{
		if(freq[i]==NULL)
		{
			freq[i] = new double[nCol];
			for(int j=0 ; j<nCol ; j++)
			{
				freq[i][j]=0.0;
			}
		}
		for(int j=0;j<nCol;j++)
		{
			myfile << freq[i][j] << " ";
		}
		myfile << endl;
		//		cout << endl;
	}
	myfile << "FIN_FREQ\n";
	myfile.close();
}
FreqKmer* FreqKmer::initFromConf(string fichier)
{
	string ligne;
	ifstream file(fichier.c_str());
	getline(file,ligne);
	FreqKmer* res;
	enum state{INIT=0, SIZE=1, FREQ=2,END=3};
	enum state s=INIT;
	double val;
	int winParam,shiftParam,j;
	bool b1_list,b2_noData;
	string pathFastaParam,pathPatternParam,keyParam,rootParam;
	while (file)
	{
		;
		if(ligne=="FIN_CONS")
		{
			s=SIZE;
		}
		if(ligne=="FIN_DIM")
		{
			s=FREQ;
		}
		if(ligne=="FIN_FREQ")
		{
			s=END;
		}
		switch (s) {
		case INIT:
			if(ligne=="1")
			{
				getline(file,ligne);
				winParam=atoi(ligne.c_str());
				getline(file,ligne);
				shiftParam=atoi(ligne.c_str());
				getline(file,ligne);
				(ligne=="0") ? (b1_list=false) : (b1_list=true);
				getline(file,ligne);
				pathFastaParam = ligne;
				getline(file,ligne);
				pathPatternParam = ligne;
				getline(file,ligne);
				(ligne=="0") ? (b2_noData=false) : (b2_noData=true);
				getline(file,ligne);
				keyParam=ligne;
				res = new FreqKmer(winParam,shiftParam,b1_list, pathFastaParam,pathPatternParam,b2_noData,keyParam);
				res->initFreq();
				getline(file,ligne);
			}
			if(ligne=="2")
			{
				getline(file,ligne);
				winParam=atoi(ligne.c_str());
				getline(file,ligne);
				(ligne=="0") ? (b1_list=false) : (b1_list=true);
				getline(file,ligne);
				pathFastaParam = ligne;
				getline(file,ligne);
				pathPatternParam = ligne;
				getline(file,ligne);
				(ligne=="0") ? (b2_noData=false) : (b2_noData=true);
				getline(file,ligne);
				keyParam=ligne;
				res = new FreqKmer(winParam,b1_list, pathFastaParam,pathPatternParam,b2_noData,keyParam);
				res->initFreq();
				getline(file,ligne);
			}
			if(ligne=="3")
			{
				getline(file,ligne);
				winParam=atoi(ligne.c_str());
				getline(file,ligne);
				pathPatternParam=ligne;
				getline(file,ligne);
				(ligne=="0") ? (b2_noData=false) : (b2_noData=true);
				getline(file,ligne);
				rootParam=ligne;
				getline(file,ligne);
				keyParam=ligne;
				res = new FreqKmer(winParam,pathPatternParam,b2_noData,rootParam,keyParam);
				res->initFreq();
				getline(file,ligne);
			}
			if(ligne=="4")
			{
				getline(file,ligne);
				winParam=atoi(ligne.c_str());
				getline(file,ligne);
				shiftParam=atoi(ligne.c_str());
				getline(file,ligne);
				pathPatternParam=ligne;
				getline(file,ligne);
				(ligne=="0") ? (b2_noData=false) : (b2_noData=true);
				getline(file,ligne);
				rootParam=ligne;
				getline(file,ligne);
				keyParam=ligne;
				res = new FreqKmer(winParam,shiftParam,pathPatternParam,b2_noData,rootParam,keyParam);
				res->initFreq();
				getline(file,ligne);
			}
			break;
		case SIZE:
			getline(file,ligne);
			res->nLine=atoi(ligne.c_str());
			getline(file,ligne);
			res->nCol=atoi(ligne.c_str());
			getline(file,ligne);
			break;
		case FREQ:
		{
			for(int i=0;i<res->nLine;i++)
			{
				if(res->freq[i]==NULL)
				{
					//cout << "la ligne " << i << " etait null \n";
					res->freq[i] = new double[res->nCol];
					for(int j=0 ; j<res->nCol ; j++)
					{
						res->freq[i][j]=0.0;
					}
				}
				getline(file,ligne);
				string buf="";
				stringstream ss(ligne);
				j=0;
				while (ss >> buf)
				{
					val=atof(buf.c_str());
					res->freq[i][j]=val;
					j+=1;
				}
			}
			s=END;
		}
		break;
		case END:
			getline(file,ligne);
			break;
		default:
			break;
		}
	}
	return res;
}
bool FreqKmer::equal(FreqKmer *f)
{
	bool res=true;
	if((f->getNCol()==this->getNCol()) && (f->getNLine()==this->getNLine()))
	{
		for(int i=0;i<f->getNLine();i++)
		{
			for(int j=0;j<f->getNCol();j++)
			{
				//                		cerr << "f[" << i << "][" << j << "] = " << f->getFreq()[i][j] <<  " || g[" << i << "][" << j << "] = " << g->getFreq()[i][j]<< "\n";
				if(f->getFreq()[i][j]!=this->getFreq()[i][j])
				{
					res=false;
				}
				if(!res)
					break;
			}
			if(!res)
				break;
		}
	}
	else
	{
		res=false;
	}
	return res;
}
void FreqKmer::writeCrossVal(int percent)
{
	writeCrossVal(percent,"learn.arff","toPredict.arff");
}
void FreqKmer::writeCrossVal(int percent,string outLearn,string outToclassify)
{
	if (dataVerbose)
	{
		cerr << "Debut FreqKmer::writeCrossVal("<< percent << ")\n";
	}
	int nbSeq = getNbAllTrue();
	int nbToTake;
	int cpt=-1;
	int seqInData_i;
	vector<int> candidates;
	nbToTake = (percent*nbSeq)/100;
	if (pathRoot!="null")
	{
		outLearn = pathRoot+"/frequencies/"+outLearn;
		outToclassify = pathRoot+"/frequencies/"+outToclassify;
	}
	ofstream os_learn ;
	ofstream os_predict ;
	bool *seqInLearn=NULL;
	if (nbToTake==0)
	{
		cerr << "WARNING in FreqKmer::writeCrossVal, percent =  " << percent << " is too low (total seq = " << nbSeq << ")\n";
		exit(0);
	}
	else
	{
		os_learn.open(outLearn.c_str());
		os_predict.open(outToclassify.c_str());
		writeHeaderWeka(os_predict);
		writeHeaderWeka(os_learn);
		seqInLearn = new bool[nbSeq];
		for(int i=0;i<nbSeq;i++)
		{
			seqInLearn[i]=true;
		}
		randomTab(&candidates,nbSeq,nbToTake);
		for(unsigned j=0;j<candidates.size();j++)
		{
			seqInLearn[candidates[j]]=false;
		}
	}
	// Pour chaque data
	for(int i=0;i<nbFastaFile;i++)
	{

		seqInData_i= this->nbSeq[i];

		//Pour chaque seq
		for(int j=0;j<seqInData_i;j++)
		{
			//si mask[i][j]
			if(mask[i][j])
			{
				cpt++;
				//si seqInlearn[cpt]
				if(seqInLearn[cpt])
				{
					//impression dans os_learn
					writeLineInOs(os_learn,i,j);
				}
				else
				{
					//sinon impression dans predict
					writeLineInOs(os_predict,i,j);
				}
			}
		}
	}
	os_learn.close();
	os_predict.close();
	if (dataVerbose)
	{
		cerr << "Fin FreqKmer::writeCrossVal("<< percent << ")\n";
	}
}
void FreqKmer::writeHeaderWeka(ofstream &os)
{
	if (dataVerbose)
	{
		cerr << "Debut FreqKmer::writeHeaderWeka(ofstream &os)\n";
	}
	os << "@RELATION freqKmer\n\n";
	vector<string> combi;
	for(int i=0;i<getNPattern();i++)
	{
		combi = patterns[i]->getCombi();
		for(unsigned j=0;j<combi.size();j++)
		{
			os << "@ATTRIBUTE " << combi[j] << " NUMERIC\n";
		}
	}
	os << "@ATTRIBUTE class {";
	for(int i=0;i<nbChildTaxa-1;i++)
	{
		os << getIdTaxa(i) << ",";
	}
	os << getIdTaxa(nbChildTaxa-1) << "}\n\n\n@DATA\n";
	if (dataVerbose)
	{
		cerr << "Fin FreqKmer::writeHeaderWeka(ofstream &os)\n";
	}
}
void FreqKmer::writeLineInOs(ofstream &os,int i,int j)
{
	if (dataVerbose)
	{
		cerr << "Debut FreqKmer::writeLineInOs\n";
	}
	int start = obtainStartLineDataSeq(i,j);
	int end = obtainEndLineDataSeq(i,j);
	string taxid = idTaxaFromData[i];
	//if(bigData)
	if(!freqFilled)
	{
		initFreq();
		fillFreq(i,j);
	}
	for(int l=start; l <= end ;l++)
	{
		//		os << "(" << i << "," << j << "):" ;
		if(normalizeBool)
		{
			for(int t=0;t<nPattern;t++)
			{
				normalizeLine(t,l);
			}
		}

		if(freq[l]!=NULL)
		{
			for(int k=0;k<nCol;k++)
			{
				os << freq[l][k] << ",";
			}
			os << taxid << " % Data(" << i << "," << j << ")\n";
		}
	}
	if(!freqFilled)
	{
		for (int var=start; var <= end ; var++)
		{
			if(freq[var]!=NULL)
			{
				//			cout << "DELETE freq[" << var << "]\n";
				delete[] freq[var];
				freq[var]=NULL;
			}
		}
	}
	if (dataVerbose)
	{
		cerr << "Fin FreqKmer::writeLineInOs\n";
	}
}
string FreqKmer::getTaxidFromString(string line)
{
	string res="";
	int index ;
	string tmp="";
	index = line.find("/data/fasta/");
	tmp=line.substr(0,index);
	index = tmp.find_last_of("__");
	res = tmp.substr(index+1,tmp.size());
	/* on chope le other du taxon courant taxon_courant/others si il y a lieu*/
	if (res.find("others") != std::string::npos)
	{
		res="others";
	}
	return res;
}
void FreqKmer::fillFreq(int data_i,int seq_j)
{
	if(dataVerbose){
		cerr << "FreqKmer::fillFreq(int data_i,int seq_j)\n ";
		cerr.flush();
	}
	int start_line = -1;
	int end_line;
	if(!freqFilled)
	{
		/* init les cases à 0 */
		/* Pour chaque kmer */
		for(int k=0;k<nPattern;k++)
		{
			/* numéro de ligne */
			if(dataVerbose)
			{
				cerr << "--------Traitement Data["<< data_i << "]\n ";
				cerr.flush();
			}
			if(noData)
			{
				if(dataVerbose)
				{
					cerr << "----------------noData=true new Data()\n ";
					cerr.flush();
				}
				if(data[data_i]!=NULL)
				{
					cerr << "WARNING in FreqKmer::fillFreq(), NULL pointer on data[" << data_i << "] expected\n";
					exit(0);
				}
				else
				{
					data[data_i] = new Data();
					//					cerr << "debut init de " << pathFasta[i] << "\n";
					data[data_i]->initFrom(pathFasta[data_i],Fasta);
					//					cerr << "fin init \n";
				}
			}
			if(dataVerbose)
			{
				cerr << "----------------Traitement Data["<< data_i << "]["<< seq_j << "]\n ";
				cerr.flush();
			}
			/*count(int *seq,int tailleDeLaSequence,int indiceKmer) */
			if(takeDataSeq(data_i,seq_j))
			{
				start_line = obtainStartLineDataSeq(data_i,seq_j);
				end_line = obtainEndLineDataSeq(data_i,seq_j);
				for(int w=start_line;w<=end_line;w++)
				{
					if(freq[w]==NULL)
					{
						//						cout << "ALLOC freq[" << w << "]\n";
						freq[w] = new double[nCol];
						for(int x=0;x<nCol;x++)
						{
							freq[w][x]=0.0;
						}
					}
				}
				if(data[data_i]==NULL)
				{
					cerr << "WARNING in FreqKmer::fillFreq(int data_i,int seq_j), NULL pointer on data[" << data_i << "]\n";
					exit(0);
				}
				count(data[data_i]->getDataObject()[seq_j],data[data_i]->getLengthSeq(seq_j),k,start_line);
			}
			if(noData)
			{
				delete data[data_i];
				data[data_i]=NULL;
			}
		}
	}
	if(dataVerbose){
		cerr << "FIN  FreqKmer::fillFreq(int data_i,int seq_j)\n ";
		cerr.flush();
	}
}
void FreqKmer::writeCrossVal(int percent,int id)
{
	stringstream sstm1,sstm2;
	sstm1 << "learn-" << id << ".arff";
	sstm2 << "toPredict-" << id << ".arff";
	string learn = sstm1.str();
	string toPredict = sstm2.str();
	writeCrossVal(percent,learn,toPredict);
}
FreqKmer* FreqKmer::sampleMe(vector<pair<int, int> > list)
{
	FreqKmer *res;
	if (dataVerbose)
	{
		cerr << "Debut FreqKmer::sampleMe(listTaxon)\n";
	}
	if(initFromRoot)
	{
		if(initWithJump)
		{
			res = new FreqKmer(winSize,shift,pathPattern,noData,pathRoot,key_fasta);
		}
		else
		{
			res = new FreqKmer(winSize,pathPattern,noData,pathRoot,key_fasta);
		}
	}
	else
	{
		if(initWithJump)
		{
			res = new FreqKmer(winSize,shift,isList,pathFastaFile,pathPattern,noData,key_fasta);
		}
		else
		{
			res = new FreqKmer(winSize,isList,pathFastaFile,pathPattern,noData,key_fasta);
		}
	}

	if(LOADALL && noData)
	{
		res->setNoData(false);
		res->loadAll();
	}

	res->freqFilled=this->freqFilled;
	/* on recupére les lignes calculées */
	if(this->freqFilled)
	{
		res->initFreq();
		for(int i=0;i<nLine;i++)
		{
			if(res->freq[i]==NULL)
			{
				res->freq[i] = new double[nCol];
				for(int p=0 ; p<res->nCol ; p++)
				{
					res->freq[i][p]=0.0;
				}
			}
			for(int j=0;j<nCol;j++)
			{
				res->freq[i][j]=this->freq[i][j];
			}
		}
	}
	int nbSequences;
	for(int i=0;i<nbFastaFile;i++)
	{

		nbSequences = nbSeq[i];
		for(int j=0;j<nbSequences;j++)
		{
			res->mask[i][j]=false;
		}
	}
	/* MAJ de la table mask */
	for(unsigned int l=0;l<list.size();l++)
	{
		res->mask[list[l].first][list[l].second]=true;
		//		cout << "mask[" << list[l].first << "][" << list[l].second << "] = true \n";
	}
	if (dataVerbose)
	{
		cerr << "Fin FreqKmer::sampleMe(list)\n";
	}
	return res;
}
void FreqKmer::writeCrossVal(FreqKmer *freqLearn, FreqKmer *freqPredict, int percent, string outLearn, string outToclassify)
{
	if (dataVerbose)
	{
		cerr << "Debut FreqKmer::writeCrossVal("<< percent << ")\n";
	}
	int nbSeq = getNbAllTrue();
	int nbToTake;
	int cpt=-1;
	int seqInData_i;
	vector<int> candidates;
	nbToTake = (percent*nbSeq)/100;
	if (pathRoot!="null")
	{
		outLearn = pathRoot+"/frequencies/"+outLearn;
		outToclassify = pathRoot+"/frequencies/"+outToclassify;
	}
	ofstream os_learn ;
	ofstream os_predict ;
	bool *seqInLearn=NULL;
	if (nbToTake==0)
	{
		cerr << "WARNING in FreqKmer::writeCrossVal, percent =  " << percent << " is too low (total seq = " << nbSeq << ")\n";
		exit(0);
	}
	else
	{
		os_learn.open(outLearn.c_str());
		os_predict.open(outToclassify.c_str());
		writeHeaderWeka(os_predict);
		writeHeaderWeka(os_learn);
		seqInLearn = new bool[nbSeq];
		for(int i=0;i<nbSeq;i++)
		{
			seqInLearn[i]=true;
		}
		randomTab(&candidates,nbSeq,nbToTake);
		for(unsigned j=0;j<candidates.size();j++)
		{
			seqInLearn[candidates[j]]=false;
		}
	}
	// Pour chaque data
	for(int i=0;i<nbFastaFile;i++)
	{

		seqInData_i= this->nbSeq[i];

		//Pour chaque seq
		for(int j=0;j<seqInData_i;j++)
		{
			//si mask[i][j]
			if(mask[i][j])
			{
				cpt++;
				//si seqInlearn[cpt]
				if(seqInLearn[cpt])
				{
					//impression dans os_learn
					writeLineInOs(os_learn,i,j,freqLearn);
				}
				else
				{
					//sinon impression dans predict
					writeLineInOs(os_predict,i,j,freqPredict);
				}
			}
		}
	}
	os_learn.close();
	os_predict.close();
	if (dataVerbose)
	{
		cerr << "Fin FreqKmer::writeCrossVal("<< percent << ")\n";
	}
	if(seqInLearn!=NULL)
		delete[] seqInLearn;
}
void FreqKmer::writeCrossVal(FreqKmer *freqPredict, int percent, string outToclassify)
{
	if (dataVerbose)
	{
		cerr << "Debut FreqKmer::writeCrossVal("<< percent << ")\n";
	}
	int nbSeq = getNbAllTrue();
	int nbToTake;
	int cpt=-1;
	int seqInData_i;
	vector<int> candidates;
	nbToTake = (percent*nbSeq)/100;
	if (pathRoot!="null")
	{
		outToclassify = pathRoot+"/frequencies/"+outToclassify;
	}
	ofstream os_predict ;
	bool *seqInLearn=NULL;
	if (nbToTake==0)
	{
		cerr << "WARNING in FreqKmer::writeCrossVal, percent =  " << percent << " is too low (total seq = " << nbSeq << ")\n";
		exit(0);
	}
	else
	{
		os_predict.open(outToclassify.c_str());
		writeHeaderWeka(os_predict);
		seqInLearn = new bool[nbSeq];
		for(int i=0;i<nbSeq;i++)
		{
			seqInLearn[i]=true;
		}
		randomTab(&candidates,nbSeq,nbToTake);
		for(unsigned j=0;j<candidates.size();j++)
		{
			seqInLearn[candidates[j]]=false;
		}
	}
	// Pour chaque data
	for(int i=0;i<nbFastaFile;i++)
	{

		seqInData_i= this->nbSeq[i];

		//Pour chaque seq
		for(int j=0;j<seqInData_i;j++)
		{
			//si mask[i][j]
			if(mask[i][j])
			{
				cpt++;
				//si seqInlearn[cpt]
				if(!seqInLearn[cpt])
				{
					//impression dans predic_t
					writeLineInOs(os_predict,i,j,freqPredict);
				}
			}
		}
	}
	os_predict.close();
	if (dataVerbose)
	{
		cerr << "Fin FreqKmer::writeCrossVal("<< percent << ")\n";
	}
	if(seqInLearn!=NULL)
		delete[] seqInLearn;
}
void FreqKmer::writeLineInOs(ofstream &os,int i,int j,FreqKmer *f)
{
	if (dataVerbose)
	{
		cerr << "Debut FreqKmer::writeLineInOs(ofstream &os,int i,int j,FreqKmer *f)\n";
	}
	int start = f->obtainStartLineDataSeq(i,j);
	int end = f->obtainEndLineDataSeq(i,j);
	string taxid = f->idTaxaFromData[i];
	//if(bigData)
	if(!f->freqFilled)
	{
		if(f->freq==NULL)
			f->initFreq();
		f->fillFreq(i,j);
	}
	for(int l=start; l <= end ;l++)
	{
		//		os << "(" << i << "," << j << "):" ;
		if(normalizeBool)
		{
			for(int t=0;t<f->nPattern;t++)
			{
				f->normalizeLine(t,l);
			}
		}
		if(f->freq[l]!=NULL)
		{
			for(int k=0;k<f->nCol;k++)
			{
				os << f->freq[l][k] << ",";
			}
			os << taxid << " % Data(" << i << "," << j << ")\n";
		}
	}
	if(!f->freqFilled)
	{
		for (int var=start; var <= end ; var++)
		{
			if(f->freq[var]!=NULL)
			{
				//			cout << "DELETE freq[" << var << "]\n";
				delete[] f->freq[var];
				f->freq[var]=NULL;
			}
		}
	}
	if (dataVerbose)
	{
		cerr << "Fin FreqKmer::writeLineInOs(ofstream &os,int i,int j,FreqKmer *f)\n";
	}
}
void FreqKmer::writeNCrossVal(FreqKmer *freqLearn, FreqKmer *freqPredict, int percent, int i,string id)
{
	stringstream sstm1,sstm2;
	int size = freqLearn->winSize;
	string s;
	if(size==-1)
	{
		s="complete";
	}
	else
	{
		s = to_string(size);
	}
	sstm1 << s << "_learn-" << i << ".arff";
	sstm2 << id << "_toPredict-" << i << ".arff";
	string learn = sstm1.str();
	string toPredict = sstm2.str();
	writeCrossVal(freqLearn,freqPredict,percent,learn,toPredict);
}
void FreqKmer::writeNCrossVal(FreqKmer *freqLearn, FreqKmer *freqPredict, int percent, int i)
{
	writeNCrossVal(freqLearn,freqPredict,percent,i,"");
}
void FreqKmer::writeNCrossVal(FreqKmer *freqPredict, int percent, int i,string id)
{
	stringstream sstm2;
	sstm2 << id << "_toPredict-" << i << ".arff";
	string toPredict = sstm2.str();
	writeCrossVal(freqPredict,percent,toPredict);
}
void FreqKmer::generateWekaData(int sizeSample,int percent,int start_win_predict,int end, int pas, int nCross)
{
	FreqKmer *toLearn = NULL;
	FreqKmer *toPredict = NULL;
	FreqKmer *res ;
	std::string s = to_string(start_win_predict);

	if(start_win_predict==-1)
	{
		s="complete";
	}
	int size = start_win_predict;
	if(initFromRoot)
	{
		cout << "new res\n";
		if(initWithJump)
		{
			res = new FreqKmer(start_win_predict,shift,pathPattern,noData,pathRoot,key_fasta);
		}
		else
		{
			res = new FreqKmer(start_win_predict,pathPattern,noData,pathRoot,key_fasta);
		}
	}
	else
	{
		if(initWithJump)
		{
			res = new FreqKmer(start_win_predict,shift,isList,pathFastaFile,pathPattern,noData,key_fasta);
		}
		else
		{
			res = new FreqKmer(start_win_predict,isList,pathFastaFile,pathPattern,noData,key_fasta);
		}
	}
	res->freqFilled=this->freqFilled;
	/* on recupére les lignes calculées */
	if(this->freqFilled)
	{
		res->initFreq();
		for(int i=0;i<nLine;i++)
		{
			if(res->freq[i]==NULL)
			{
				res->freq[i] = new double[nCol];
				for(int p=0 ; p<res->nCol ; p++)
				{
					res->freq[i][p]=0.0;
				}
			}
			for(int j=0;j<nCol;j++)
			{
				res->freq[i][j]=this->freq[i][j];
			}
		}
	}
	if(sizeSample>0)
	{
		toLearn = sampleMe(sizeSample);
		toPredict = res->sampleMe(getSampledTaxon());
	}
	else
	{
		toLearn = sampleMe(-1);
		toPredict = sampleMe(-1);
	}
	for(int i=1;i<=nCross;i++)
	{
		cout << "writeNCrossVal(" << i << ")\n";
		writeNCrossVal(toLearn,toPredict,percent,i,s);
	}
	size += pas;
	cout << "delete res\n";
	delete res;
	delete toPredict;
	cout << "toto\n";
	delete toLearn;


	if(start_win_predict>0)
	{
		for(int i=size; i<=end;i+=pas)
		{
			s = to_string(i);
			cout << "new res tailleF = " << s << "\n";
			if(initFromRoot)
			{
				if(initWithJump)
				{
					res = new FreqKmer(i,shift,pathPattern,noData,pathRoot,key_fasta);
				}
				else
				{
					res = new FreqKmer(i,pathPattern,noData,pathRoot,key_fasta);
				}
			}
			else
			{
				if(initWithJump)
				{
					res = new FreqKmer(i,shift,isList,pathFastaFile,pathPattern,noData,key_fasta);
				}
				else
				{
					res = new FreqKmer(i,isList,pathFastaFile,pathPattern,noData,key_fasta);
				}
			}
			if(sizeSample>0)
			{
				toPredict = res->sampleMe(getSampledTaxon());
			}
			else
			{
				toPredict = sampleMe(-1);
			}
			for(int j=1;j<=nCross;j++)
			{
				cout << "writeNCrossVal(" << j << ")\n";
				writeNCrossVal(toPredict,percent,j,s);
			}
			cout << "delete res\n";
			delete res;
			delete toPredict;
		}
	}
}

void FreqKmer::setNoData(bool b)
{
	noData=b;
}

int FreqKmer::getNbSeqData(int data_i)
{
	return nbSeq[data_i];
}

bool FreqKmer::existSeqTrueInData(int data_i)
{
	bool res = false;
	for(int i=0;i<getNbSeqData(data_i);i++)
	{
		if(takeDataSeq(data_i,i))
		{
			res = true;
		}
		if (res)
			break;
	}

	return res;
}
void FreqKmer::loadAll()
{
	if(noData==false)
	{
		for(int i=0;i<nbFastaFile;i++)
		{
			if(existSeqTrueInData(i))
			{
				if(data[i]!=NULL)
				{
					cerr << "WARNING in FreqKmer::loadAll(), NULL pointer on data[" << i << "] expected\n";
					exit(0);
				}
				data[i] = new Data();
				data[i]->initFrom(pathFasta[i],Fasta);
			}
		}
	}
}
