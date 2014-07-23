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

FreqKmer::FreqKmer() {
	/* On veut les accessions donc Yes */
	data=NULL;
	freq=NULL;
	patterns=NULL;
	kmerSpace=NULL;
	indexLineDataSeq=NULL;
	indexLineData=NULL;
	nCol=0;
	nSeq=0;
	nLine=0;
	nPattern=0;
	nbFichierFasta=0;
	/* Taille par défaut */
	winSize=WINDOW_SIZE;
	index=0;
	shift=1;
}

FreqKmer::FreqKmer(int win_size, int s) {
	data=NULL;
	freq=NULL;
	patterns=NULL;
	kmerSpace=NULL;
	shift=s;
	nCol=0;
	nSeq=0;
	nLine=0;
	nPattern=0;
	nbFichierFasta=0;
	winSize=win_size;
	index=0;
	indexLineDataSeq=NULL;
	indexLineData=NULL;

}

FreqKmer::FreqKmer(int win_size) {
	data=NULL;
	freq=NULL;
	patterns=NULL;
	kmerSpace=NULL;
	nCol=0;
	nSeq=0;
	nLine=0;
	nPattern=0;
	nbFichierFasta=0;
	winSize=win_size;
	index=0;
	shift=1;
	indexLineDataSeq=NULL;
	indexLineData=NULL;

}

FreqKmer::~FreqKmer() {

	if(dataVerbose){
		cerr << "Debut FreqKmer::~FreqKmer()\n ";
		cerr.flush();
	}
	for (int var = 0; var < nPattern ; var++)
	{
		delete patterns[var];
	}
	delete[] patterns;
	for (int var=0; var < nbFichierFasta ; var++)
	{
		delete data[var];
		delete indexLineDataSeq[var];
	}
	delete[] data;
	delete[] indexLineDataSeq;

	if(dataVerbose){
		cerr << "Debut delete freq\n ";
		cerr.flush();
	}
	if (freq!=NULL)
	{

		for (int var=0; var < nLine ; var++)
		{
			delete freq[var];
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
	while (file)
	{
		if(tailleLigne!=0)
		{
			nbFichier++;
		}
		getline(file,ligne);
	}

	nbFichierFasta=nbFichier;
	data=new Data*[nbFichier];
	indexLineData = new int[nbFichier];
	ifstream file2(fichier.c_str());
	getline(file2,ligne);
	/* J'init à partir de chaque ligne du fichier */
	indexLineDataSeq = new int*[nbFichierFasta];
	while (file2)
	{
		if(tailleLigne!=0)
		{
			data[cpt] = new Data(Yes);
			//			cerr << "lecture de " << ligne << "\n";
			data[cpt]->initFrom(ligne,Fasta);
			/* Ajouter des attributs ??*/
			indexLineDataSeq[cpt] = new int[data[cpt]->getNtaxa()];
			nSeq += data[cpt]->getNtaxa();
			if (winSize>0)
			{
				for (int var = 0; var < data[cpt]->getNtaxa(); var++)
				{

					tailleSeq = data[cpt]->getPrimarySequence(var).length();

					/* Fenetre plus grande que la sequence */
					if (tailleSeq<winSize)
					{
						nLine += 1;
					}
					else
					{
						//nLine += tailleSeq-winSize+1;
						nLine += obtainNbLineWindow(0,tailleSeq-1,winSize,shift);
					}

					//					cerr << "indexLineDataSeq[cpt][var] <-- " << nLine << "\n";
					indexLineDataSeq[cpt][var]=nLine-1;
				}
			}
			else
			{
				nLine +=  data[cpt]->getNtaxa();
				for(int i=0;i<data[cpt]->getNtaxa();i++)
				{
					indexLineDataSeq[cpt][i] = cpt2++;
				}

			}

			indexLineData[cpt]=nLine-1;
			cpt++;
		}
		getline(file2,ligne);
	}
}

void FreqKmer::initFromFasta(string fichier)
{
	nbFichierFasta=1;
	data=new Data*[1];
	indexLineData = new int[1];
	int cpt2=-1;

	data[0] = new Data(Yes);
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
}

int FreqKmer::obtainColIndex(int indexPattern,int *seq,int pos)
{
	//	int tmp = patterns[indexPattern]->getKmer(seq,pos);
	//	int s=obtainStartColKmer(indexPattern);
	//	return s+tmp;

	return obtainStartColKmer(indexPattern)+patterns[indexPattern]->getKmer(seq,pos);
}


void FreqKmer::initFreq()
{
	if(dataVerbose){
		cerr << "Debut FreqKmer::initFreq()\n ";
		cerr.flush();
	}
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
	int col;
	int cpt=0;
	/* on s'arrête à la dernière fenetre possible */
	int fin = pos+win_length-patterns[indexPattern]->getTaillePattern();
	for(int i = pos; i <= fin ; i++)
	{
		/* indice du kmer à la position i */
		col = obtainColIndex(indexPattern,seq,i);
		freq[index][col]+=1;
		previous[cpt]=col;
		//		cerr << "previous[" << i << "] reçoit " << col <<"\n";
		cpt++;
	}


}

void FreqKmer::swapBuffAndCount(int *current,int *previous,int buf_size, int indexPattern,int *seq,int pos)
{
	for(int i=shift;i<buf_size;i++)
	{
		current[i-shift]=previous[i];
	}
	for(int i=buf_size;i<buf_size+shift;i++)
	{
		current[i-shift] = obtainColIndex(indexPattern,seq,pos+i-shift);
	}
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

	/*taille de la sous sequence où travailler */
	int taille_sous_sequence = 0;
	int z;

	/* Initialisation
	 * des tableaux premier
	 * et dernier
	 */


	if (winSize==-1 || winSize>seq_length)
	{
		taille_sous_sequence = seq_length;
	}
	else
	{
		taille_sous_sequence = winSize;
	}

	int buf_size = taille_sous_sequence-patterns[indexPattern]->getTaillePattern()+1;
	int *previous = new int[buf_size];
	int *current = new int[buf_size];

	/* on recupère le nombre de ligne pour la sequence
	 * afin d'iterer le bon nombre de fois
	 */
	z=obtainNbLineWindow(0,seq_length-1,taille_sous_sequence,shift);

	while(j < z)
	{

		/* winCount(seq,winSize,indiceDepart,indice pattern */
		if (j>0)
		{


			swapBuffAndCount(current,previous,buf_size,indexPattern,seq,i);
			for(int i=0;i<buf_size;i++)
			{
				freq[index][current[i]]+=1;
			}
			swap(current,previous,buf_size);



		}
		else
		{
			winCount(seq,taille_sous_sequence,i,indexPattern,previous);

		}



		/* On decale de shift nucle */

		i = i + shift;
		j = j + 1;
		index++;
	}
	//
	//	cerr << "      i = " << i << "\n";
	//	cerr << " nLine = " << nLine << "\n";
	//	cerr << " index  = " << index << "\n";
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
		for(int i=0;i<nbFichierFasta;i++)
		{
			if(dataVerbose)
			{
				cerr << "--------Traitement Data["<< i << "]\n ";
				cerr.flush();
			}
			/* Pour chaque Data du fichier */
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
			//			cout << freq[i][j] << ";";
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
	//	int res = 0;
	//
	//	do
	//	{
	//		res += 1;
	//		i += pas;
	//
	//	} while (i<=(j-l+1));
	//
	//	return res;
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
	if(j==0)
	{
		if(i==0)
		{
			return 0;
		}
		else
		{
			return indexLineDataSeq[i-1][data[i-1]->getNtaxa()-1]+1;
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
