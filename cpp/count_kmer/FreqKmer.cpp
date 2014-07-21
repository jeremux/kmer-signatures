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
	nData=0;
	nLigne=0;
	nPattern=0;
	nbFichierFasta=0;
	/* Taille par défaut */
	tailleFenetre=TAILLE_FENETRE;
	index=0;
	shift=1;
	tabDernier=new int[shift];
	tabPremier=new int[shift];
}

FreqKmer::FreqKmer(int tailleF, int s) {
	data=NULL;
	freq=NULL;
	patterns=NULL;
	kmerSpace=NULL;
	shift=s;
	nCol=0;
	nData=0;
	nLigne=0;
	nPattern=0;
	nbFichierFasta=0;
	tailleFenetre=tailleF;
	index=0;
	tabDernier=new int[shift];
	tabPremier=new int[shift];
	indexLineDataSeq=NULL;
	indexLineData=NULL;

}

FreqKmer::FreqKmer(int tailleF) {
	data=NULL;
	freq=NULL;
	patterns=NULL;
	kmerSpace=NULL;
	nCol=0;
	nData=0;
	nLigne=0;
	nPattern=0;
	nbFichierFasta=0;
	tailleFenetre=tailleF;
	index=0;
	shift=1;
	tabPremier=new int[shift];
	tabDernier=new int[shift];
	indexLineDataSeq=NULL;
	indexLineData=NULL;

}

FreqKmer::~FreqKmer() {

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

	if (freq!=NULL)
	{
		for (int var=0; var < nLigne ; var++)
		{
			delete freq[var];
		}
		delete[] freq;
	}
	delete[] kmerSpace;
	delete[] tabPremier;
	delete[] tabDernier;
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
			nData += data[cpt]->getNtaxa();
			if (tailleFenetre>0)
			{
				for (int var = 0; var < data[cpt]->getNtaxa(); var++)
				{

					tailleSeq = data[cpt]->getPrimarySequence(var).length();

					/* Fenetre plus grande que la sequence */
					if (tailleSeq<tailleFenetre)
					{
						nLigne += 1;
					}
					else
					{
						//nLigne += tailleSeq-tailleFenetre+1;
						nLigne += obtainNbLineWindow(0,tailleSeq-1,tailleFenetre,shift);
					}

//					cerr << "indexLineDataSeq[cpt][var] <-- " << nLigne << "\n";
					indexLineDataSeq[cpt][var]=nLigne-1;
				}
			}
			else
			{
				nLigne +=  data[cpt]->getNtaxa();
				for(int i=0;i<data[cpt]->getNtaxa();i++)
				{
					indexLineDataSeq[cpt][i] = cpt2++;
				}

			}

			indexLineData[cpt]=nLigne-1;
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

	nData = data[0]->getNtaxa();
	indexLineDataSeq = new int *[1];

	indexLineDataSeq[0] = new int[nData];

	int taille=0;
	if (tailleFenetre>0)
	{
		for (int var = 0; var < nData; var++)
		{
			taille= data[0]->getPrimarySequence(var).length();
			if (taille<tailleFenetre)
			{
				nLigne += 1;
			}
			else
			{
				//nLigne += taille-tailleFenetre+1;
				nLigne += obtainNbLineWindow(0,taille-1,tailleFenetre,shift);
			}
			indexLineDataSeq[0][var]=nLigne-1;
		}
	}
	else
	{
		for (int var = 0; var < nData; var++)
		{
			indexLineDataSeq[0][var]=cpt2++;
		}
		nLigne = nData;
	}
	indexLineData[0] = nLigne-1;
}

int FreqKmer::obtainColIndex(int indicePattern,int *seq,int pos)
{
//	int tmp = patterns[indicePattern]->getKmer(seq,pos);
//	int s=obtainStartColKmer(indicePattern);
//	return s+tmp;

	return obtainStartColKmer(indicePattern)+patterns[indicePattern]->getKmer(seq,pos);
}
/**
 * Permet de copier une ligne
 * et mettant a jour les frequences
 * dans le cas d'une graine continue
 */
void FreqKmer::copieLigneFreq(int src,int dest,int indicePattern)
{
	int decalage = 0;

	for(int i=0 ; i < indicePattern ; i++)
	{
		decalage += patterns[i]->getAllCombi();
	}

//	cerr << "on va copier de " << decalage << " a " << patterns[indicePattern]->getAllCombi()+decalage-1 << "\n";
	for(int i=decalage; i<patterns[indicePattern]->getAllCombi()+decalage ; i++)
	{
		freq[dest][i]=freq[src][i];
	}

	decrementPremier(dest);
}

void FreqKmer::decrementPremier(int ligne)
{
	for(int i=0;i<shift;i++)
	{

		freq[ligne][tabPremier[i]] = freq[ligne][tabPremier[i]] - 1;
	}
}

void FreqKmer::initFreq()
{
	freq = new double *[nLigne];
	for(int i=0 ; i<nLigne ; i++)
	{
		freq[i] = new double[nCol];
		for(int j=0 ; j<nCol ; j++)
		{
			freq[i][j]=0.0;
		}
	}
}

/***
 * effectue le comptage dans une sous sequence
 * @param	seq: sequence où compter
 * @param	seq_taille: taille de la sous sequence
 * @param	debut: où commencer à compter
 * @param	indicePattern: indice du kmer courant
 */
void FreqKmer::compteFenetre(int *seq,int seq_taille,int debut ,int indicePattern)
{
	int col;
	/* on s'arrête à la dernière fenetre possible */
	int fin = debut+seq_taille-patterns[indicePattern]->getTaillePattern();
	for(int i = debut; i <= fin ; i++)
	{
		/* indice du kmer à la position i */
		col = obtainColIndex(indicePattern,seq,i);

		/* si on traite les shift premiers kmer */
		if(i-debut<shift)
		{
			/* alors on les mémorise (pour ensuite décrémenter, lors du prochain decalage */
			tabPremier[i-debut] = col;
		}
		/* on recupere le dernier kmer, pour ensuite calculer les shift-1 derniers kmers restants */
		if(i==fin)
		{
			tabDernier[0] = col;
		}
		/* dans tous les cas on augmente la frequence du kmer rencontré */
		freq[index][col]+=1;
	}
}

/**
 * Permet de passer de mettre a jour la
 * table des frequence en incrementer les
 * nouveaux derniers et en mettant à jour
 * les nouveaux premier
 * @param	seq	sequence où l'on compte
 * @param	i	endroit où on se trouve dans la seq
 * @param	seq_taille	taille de la sequence (ici la fenetre)
 * @param 	indicePattern indice du pattern courant
 */
void FreqKmer::add_one(int *seq,int i,int seq_taille,int indicePattern)
{
	int pattern_taille = patterns[indicePattern]->getTaillePattern();
	int alpha =  i + seq_taille - (pattern_taille + 1);
	int beta = i + seq_taille - 1 ;
	int kmer_taille = patterns[indicePattern]->getTailleKmer();

	int decalage = 0;

	for(int j=0 ; j < indicePattern ; j++)
	{
		decalage += patterns[j]->getAllCombi();
	}

	// printf("seq[%d] = %d\n",beta,seq[beta]);
	// printf("dernier before = %d\n",dernier);
	setDernier(seq,alpha,beta,kmer_taille,indicePattern,i,pattern_taille,seq_taille,decalage);

	incrementDernier(index);


	initPremier();
	setPremier(seq,pattern_taille,indicePattern,kmer_taille,decalage,i);

}

void FreqKmer::incrementDernier(int ligne)
{
	for(int i=0;i<shift;i++)
		freq[index][tabDernier[i]] = freq[index][tabDernier[i]] + 1;
}

void FreqKmer::setDernier(int *seq,int alpha,int beta,int kmer_taille,int indicePattern,int i,int pattern_taille,int seq_taille,int d)
{

//	if (patterns[indicePattern]->isContinue())
//	{
//		for(int j=0;j<shift;j++)
//		{
//			tabDernier[j] = (tabDernier[j]-d - (seq[alpha+j]*pow(4,kmer_taille-1)))*4 + seq[beta+j];
//			tabDernier[j] += d;
//		}
//		//cerr << "dernier = " << dernier << "\n";
//	}
//	else
	{
		for(int j=0;j<shift;j++)
			tabDernier[j] = obtainColIndex(indicePattern,seq,i+seq_taille-pattern_taille-j);
	}
//	for(int j=0;j< shift ;j++)
//	{
//		cerr << "tabDernier[" << j << "] = " << tabDernier[j] << "\n";
//	}
}

void FreqKmer::setPremier(int *seq,int pattern_taille,int indicePattern,int kmer_taille,int decalage,int i)
{
	int k,l;
	for(int j=0;j<shift;j++)
	{
		for (k = 0,l=0; k < pattern_taille; k++)
		{
			if(patterns[indicePattern]->extraire(k))
			{
				tabPremier[j] += seq[i+k+j] * pow(4,(kmer_taille-l-1));
				l++;
			}
		}
		tabPremier[j] += decalage;
	}

//	for(int j=0;j< shift ;j++)
//		{
//			cerr << "tabPremier[" << j << "] = " << tabPremier[j] << "\n";
//		}
}
/**
 * Effectue le comptaga pour seq
 * @param 	seq: la sequence ou compter
 * @param 	seq_taille : taille de la sequence
 * @parm	indicePattern : indice du kmer actuel
 */
void FreqKmer::count(int *seq,int seq_taille,int indicePattern)
{
	int i=0;
	int j=0;

	/*taille de la sous sequence où travailler */
	int taille_sous_sequence = 0;
	int z;

	/* Initialisation
	 * des tableaux premier
	 * et dernier
	 */
	initPremier();
	initDernier();

	if (tailleFenetre==-1 || tailleFenetre>seq_taille)
	{
		taille_sous_sequence = seq_taille;
	}
	else
	{
		taille_sous_sequence = tailleFenetre;
	}
	
	/* on recupère le nombre de ligne pour la sequence
	 * afin d'iterer le bon nombre de fois
	 */
	z=obtainNbLineWindow(0,seq_taille-1,taille_sous_sequence,shift);

	while(j < z)
	{

//		if (j>0)
//		{
//
//			copieLigneFreq(index-1,index,indicePattern);
//			add_one(seq,i,taille_sous_sequence,indicePattern);
//
//		}
//		else
		{
			/* compteFenetre(seq,tailleFenetre,indiceDepart,indice pattern */
			compteFenetre(seq,taille_sous_sequence,i,indicePattern);

		}

		/* On decale de shift nucle */
		i = i + shift;
		j = j + 1;
		index++;
	}
//
//	cerr << "      i = " << i << "\n";
//	cerr << " nLigne = " << nLigne << "\n";
//	cerr << " index  = " << index << "\n";
}


void FreqKmer::fillFreq()
{
	/* init les cases à 0 */
	initFreq();

	/* Pour chaque kmer */
	for(int k=0;k<nPattern;k++)
	{
		/* numéro de ligne */
		index=0;

		/* Pour chaque fichier fasta */
		for(int i=0;i<nbFichierFasta;i++)
		{
			/* Pour chaque Data du fichier */
			for(int j=0;j<data[i]->getNtaxa();j++)
			{
				/*count(int *seq,int tailleDeLaSequence,int indiceKmer) */
				count(data[i]->getDataObject()[j],data[i]->getLengthSeq(j),k);
			}
		}

	}
}

void FreqKmer::initPremier()
{
	for(int i=0;i<shift;i++)
		tabPremier[i]=0;
}

void FreqKmer::initDernier()
{
	for(int i=0;i<shift;i++)
		tabPremier[i]=0;
}

void FreqKmer::imprimeCSV(string ouput)
{
	ofstream myfile;
	myfile.open(ouput.c_str());


	for(int j=0;j<nCol;j++)
	{
		myfile << j << ";";
	}
	myfile << endl;
	for(int i=0;i<nLigne;i++)
	{
		for(int j=0;j<nCol;j++)
		{
			myfile << freq[i][j] << ";";
		}
		myfile << endl;
	}
	myfile.close();
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

	do
	{
		res += 1;
		i += pas;

	} while (i<=(j-l+1));


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
