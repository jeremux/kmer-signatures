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

using namespace std;

FreqKmer::FreqKmer(int win_size,int s,bool list, string file,string patternFile, bool b)
{
    data=NULL;
    freq=NULL;
    patterns=NULL;
    kmerSpace=NULL;
    shift=s;
    nCol=0;
    nSeq=0;
    nLine=0;
    nPattern=0;
    nbFastaFile=0;
    winSize=win_size;
    index=0;
    indexLineDataSeq=NULL;
    indexLineData=NULL;
    listData="null";
    noData=b;

    initPatterns(patternFile);
    if (list)
    {
	initDataFromListFastaPath(file);
	initPathFasta(file);
	listData = file;
    }
    else
    {
    	initFromFasta(file);
    	pathFasta = new string[1];
    	pathFasta[0] = file;

    }


}



FreqKmer::FreqKmer(int win_size,bool list, string file,string patternFile, bool b)
{
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
    index=0;
    noData=b;
    
    if(win_size>0)
    {

	shift=(win_size*20)/100;
	if (shift==0)
	    shift=1;

    }
    indexLineDataSeq=NULL;
    indexLineData=NULL;
    initPatterns(patternFile);

    if (list)
    {

	initDataFromListFastaPath(file);
	initPathFasta(file);
	listData = file;

    }
    else
    {
	initFromFasta(file);
	pathFasta = new string[1];
	pathFasta[0] = file;
	listData = "null";
    }
}

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

    delete[] data;

    for(int i=0;i<nbFastaFile;i++)
    {
	if(indexLineDataSeq[i]!=NULL)
	{
	    delete[] indexLineDataSeq[i];
	}
    }
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
		delete[] freq[var];
	}

	delete[] freq;
    }
    if(dataVerbose){
	cerr << "fin delete freq\n ";
	cerr.flush();
    }
    delete[] kmerSpace;

    if(dataVerbose){
	cerr << "Debut FreqKmer::~FreqKmer()\n ";
	cerr.flush();
    }
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

    kmerSpace = new int[nPattern];

    for(int i=0 ; i<nPattern ; i++)
    {
	if (i==0)
	{
	    kmerSpace[i]=patterns[i]->getAllCombi()-1;
	}
	else
	{
	    kmerSpace[i]=kmerSpace[i-1]+patterns[i]->getAllCombi();
	}
	nCol += patterns[i]->getAllCombi();
    }

}

int FreqKmer::obtainStartColKmer(int i)
{
    if (i==0)
    {
	return 0;
    }
    else
    {
	return kmerSpace[i-1]+1;
    }
}

int FreqKmer::obtainEndColKmer(int i)
{
    return kmerSpace[i];
}


void FreqKmer::initDataFromListFastaPath(string fichier)
{
    string ligne;
    int cpt=0;
    int cpt2=-1;
    int tailleSeq = 0;
    ifstream file(fichier.c_str());
    getline(file,ligne);
    int nbFichier = 0;
    int tailleLigne = ligne.length();
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
    while (file2)
    {
	if(tailleLigne!=0)
	{
	    /* Yes: je sauvegarde les accessions pour plus tard, ça pourrait être utile */
	    /* au premier tour cpt=0, ensuite 1...*/

	    data[cpt] = new Data();
	    //			cerr << "lecture de " << ligne << "\n";
	    /* J'initilise ma donnée */
	    string tmpligne=ligne;
	    data[cpt]->initFrom(tmpligne,Fasta);

	    indexLineDataSeq[cpt] = new int[data[cpt]->getNtaxa()];

	    /* On incrémente le nombre de sequence total pour le comptage */
	    nSeq += data[cpt]->getNtaxa();

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
		}

	    }
	    /* Si winSize <= 0 alors la taille de la fenetre est celle de la sequence */
	    else
	    {
		nLine +=  data[cpt]->getNtaxa();
		for(int i=0;i<data[cpt]->getNtaxa();i++)
		{
		    indexLineDataSeq[cpt][i] = cpt2++;
		}

	    }

	    indexLineData[cpt]=nLine-1;

	    if(noData)
	    {
		delete data[cpt];
		data[cpt]=NULL;
	    }
	    cpt++;


	}
	getline(file2,ligne);
	tailleLigne = ligne.length();
    }
}


void FreqKmer::initPathFasta(string fichier)
{
    string ligne;
    int cpt=0;
    ifstream file(fichier.c_str());



    getline(file,ligne);

    int tailleLigne = ligne.length();

    /* On lit la liste */

    pathFasta = new string[nbFastaFile];
    while (file)
    {

	/* Si ce n'est pas une ligne blanche, vide */
	if(tailleLigne!=0)
	{
	    pathFasta[cpt]=ligne;
	}
	/* on lit la ligne suivante */
	getline(file,ligne);
	tailleLigne = ligne.length();
	cpt++;
    }
}
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

    int taille=0;
    if (winSize>0)
    {
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
	}
    }
    else
    {
	for (int var = 0; var < nSeq; var++)
	{
	    indexLineDataSeq[0][var]=cpt2++;
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

int FreqKmer::obtainColIndex(int indexPattern,int *seq,int pos)
{

    return obtainStartColKmer(indexPattern)+patterns[indexPattern]->getKmer(seq,pos);
}


void FreqKmer::initFreq()
{
    if(dataVerbose){
	cerr << "Debut FreqKmer::initFreq()\n ";
	cerr.flush();
    }
//	cerr << "nLine = " << nLine << "\n";
//	cerr << "nCol = " << nCol << "\n";
    freq = new double *[nLine];
    for(int i=0 ; i<nLine ; i++)
    {

	freq[i] = new double[nCol];
	for(int j=0 ; j<nCol ; j++)
	{
	    freq[i][j]=0.0;
	}
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
 * @param	debut: où commencer à compter
 * @param	indexPattern: indice du kmer courant
 */
void FreqKmer::winCount(int *seq,int win_length,int pos ,int indexPattern,int *previous)
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
	    freq[index][col]+=1;
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

void FreqKmer::copyBuffAndCount(int *current,int *previous,int buf_size, int indexPattern,int *seq,int pos)
{
//	cout << "start copyBuffAndCount int buf_size=" << buf_size << " int indexPattern=" << indexPattern << " int pos=" << pos << " \n";
//	cout << "previous : ";
//	cout.flush();
    for(int i=0;i<buf_size;i++)
    {
	current[i]=-1;
//	obtainColIndex(indexPattern,seq,i);
    }
//	cout << "\n";
//	cout.flush();
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
//			cerr << "currant[" << i-shift << "] = " << obtainColIndex(indexPattern,seq,i) <<"\n";
	    current[cpt++] = obtainColIndex(indexPattern,seq,i);
	}
    }
//	cout << "current :";
//	for(int i=0;i<buf_size;i++)
//	{
//	cout << current[i] << "\t";
//	}
//	cout << "\n";
//	cout.flush();
//	cout << "fin copyBuffAndCount\n";
//	cout.flush();
}

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

void FreqKmer::printBuf(int *buf,int buf_size)
{
    cerr << "[";
    for(int i=0;i<buf_size;i++)
    {
	cerr << buf[i] << "][";
    }
}

/**
 * Effectue le comptaga pour seq
 * @param 	seq: la sequence ou compter
 * @param 	seq_length : taille de la sequence
 * @parm	indexPattern : indice du kmer actuel
 */
void FreqKmer::count(int *seq,int seq_length,int indexPattern)
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

	/* Si j > 0 alors on a déjà effectué le comtpage pour la première fenetre
	 * on utilise l'astuce de decalage
	 */
	if (j>0)
	{
//			cout << "window "<< j <<"\n";
//			cout.flush();
	    /* On copie la bonne partie du buffer et on compte les nouveaux kmers */
	    copyBuffAndCount(current,previous,buf_size,indexPattern,seq,i);

	    /* on increment les kmers trouvés */
	    for(int i=0;i<buf_size;i++)
	    {
		freq[index][current[i]]+=1;
	    }
	    /* le courant devient le precedent */
	    swap(current,previous,buf_size);



	}
	else
	{
//			cout << "first window\n";
//			cout.flush();
	    winCount(seq,win_length,i,indexPattern,previous);
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


void FreqKmer::fillFreq()
{
    if(dataVerbose){
	cerr << "FreqKmer::fillFreq()\n ";
	cerr.flush();
    }

    /* init les cases à 0 */
    initFreq();

    /* Pour chaque kmer */
    for(int k=0;k<nPattern;k++)
    {
	/* numéro de ligne */
	if(dataVerbose){
	    cerr << "--Traitement pattern "<< k+1 << "\n ";
	    cerr.flush();
	}
	index=0;

	/* Pour chaque fichier fasta */
	for(int i=0;i<nbFastaFile;i++)
	{
	    if(dataVerbose)
	    {
		cerr << "--------Traitement Data["<< i << "]\n ";
		cerr.flush();
	    }
	    /* Pour chaque Data du fichier */
	    if(noData)
	    {
		if(dataVerbose)
		{
		    cerr << "----------------noData=true new Data()\n ";
		    cerr.flush();
		}
		if(data[i]!=NULL)
		{
		    cerr << "WARNING in FreqKmer::fillFreq(), NULL pointer on data[" << i << "] expected\n";
		    exit(0);
		}
		else
		{
		    data[i] = new Data();

//					cerr << "debut init de " << pathFasta[i] << "\n";
		    data[i]->initFrom(pathFasta[i],Fasta);
//					cerr << "fin init \n";
		}

	    }

	    for(int j=0;j<data[i]->getNtaxa();j++)
	    {
		if(dataVerbose)
		{
		    cerr << "----------------Traitement Data["<< i << "]["<< j << "]\n ";
		    cerr.flush();
		}
		/*count(int *seq,int tailleDeLaSequence,int indiceKmer) */
		count(data[i]->getDataObject()[j],data[i]->getLengthSeq(j),k);
	    }
	    if(noData)
	    {
		delete data[i];
		data[i]=NULL;

	    }
	}

    }

    if(dataVerbose){
	cerr << "Fin FreqKmer::fillFreq()\n ";
	cerr.flush();
    }
}


void FreqKmer::imprimeCSV(string ouput)
{
    if(dataVerbose){
	cerr << "Debut FreqKmer::imprimeCSV("<< ouput << ")\n ";
	cerr.flush();
    }

    if(dataVerbose){
	cerr << "Debut FreqKmer::imprimeCSV("<< ouput << ")\n ";
	cerr.flush();
    }

    if(dataVerbose){
	cerr << "Ouverture srteam\n ";
	cerr.flush();
    }
    ofstream myfile ;

    if(dataVerbose){
	cerr << "Ouverture fichier\n ";
	cerr.flush();
    }

    myfile.open(ouput.c_str());

    if(dataVerbose){
	cerr << "Impression en tete\n ";
	cerr.flush();
    }
    for(int j=0;j<nCol;j++)
    {
	myfile << j << ";";
	//		cout << j << ";";
    }
    myfile << endl;
    //	cout << endl;
    for(int i=0;i<nLine;i++)
    {
	for(int j=0;j<nCol;j++)
	{
	    myfile << freq[i][j] << ";";
	    cout << freq[i][j] << ";";
	}
	myfile << endl;
	//		cout << endl;
    }
    myfile.close();
    if(dataVerbose){
	cerr << "Fin FreqKmer::imprimeCSV("<< ouput << ")\n ";
	cerr.flush();
    }
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
    return ((j-i+1-l)/pas)+1;
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

int FreqKmer::obtainEndLineData(int i)
{
    return indexLineData[i];
}

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

	    if(noData)
	    {
		if(data[tmp2]!=NULL)
		{
		    cerr << "WARNING in FreqKmer::obtainStartLineDataSeq, NULL pointer on data[" << tmp2 << "] expected\n";
		    exit(0);
		}
		data[tmp2] = new Data();
		data[tmp2]->initFrom(pathFasta[tmp2],Fasta);
	    }
	    tmp = data[tmp2]->getNtaxa()-1;
	    if(noData)
	    {
		delete data[tmp2];
		data[tmp2] = NULL;
	    }

	    return indexLineDataSeq[i-1][tmp]+1;
	}
    }
    else
    {
	return indexLineDataSeq[i][j-1]+1;
    }
    return 0;
}

int FreqKmer::obtainEndLineDataSeq(int i,int j)
{
    return indexLineDataSeq[i][j];
}

int FreqKmer::obtainNbLineDataSeq(int i,int j)
{
    return obtainEndLineDataSeq(i,j)-obtainStartLineDataSeq(i,j)+1;
}

void FreqKmer::obtainDataSeqFromLine(int line,int &idData,int &idSeq)
{
    bool res = false;
    for(int i=0;i<nbFastaFile;i++)
    {
	if(noData)
	{
	    if(data[i]!=NULL)
	    {
		cerr << "WARNING in FreqKmer::obtainDataSeqFromLine, NULL pointer on data[" << i << "] expected\n";
		exit(0);
	    }
	    data[i] = new Data();
	    data[i]->initFrom(pathFasta[i],Fasta);
	}

	for(int j=0;j<data[i]->getNtaxa();j++)
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
	if(noData)
	{
	    delete data[i];
	    data[i]=NULL;
	}
	if(res)
	    break;
    }
}

void FreqKmer::obtainDataSeqWinFromLine(int line,int &idData,int &idSeq,int &idWin)
{
    bool res = false;
    for(int i=0;i<nbFastaFile;i++)
    {
	if(noData)
	{
	    if(data[i]!=NULL)
	    {
		cerr << "WARNING in FreqKmer::obtainDataSeqFromLine, NULL pointer on data[" << i << "] expected\n";
		exit(0);
	    }
	    data[i] = new Data();
	    data[i]->initFrom(pathFasta[i],Fasta);
	}
	for(int j=0;j<data[i]->getNtaxa();j++)
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
	if(noData)
	{
	    delete data[i];
	    data[i]=NULL;
	}
	if(res)
	    break;
    }
}

