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
	premier=0;
	dernier=0;
	index=0;
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
	premier=0;
	dernier=0;
	index=0;
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

int FreqKmer::getStartColKmer(int i)
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

int FreqKmer::getEndColKmer(int i)
{
	return kmerSpace[i];
}


void FreqKmer::initFromList(string fichier)
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
						nLigne += tailleSeq-tailleFenetre+1;
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
				nLigne += taille-tailleFenetre+1;
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

int FreqKmer::getCol(int indicePattern,int *seq,int pos)
{
	int tmp = patterns[indicePattern]->getKmer(seq,pos);
	int s=0;
	for(int i=0; i<indicePattern ; i++)
	{
		s += patterns[i]->getAllCombi();
	}
	return s+tmp;
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
//	cerr << "decalage = " << decalage << "\n";
//	cerr << "on va copier de " << decalage << " a " << patterns[indicePattern]->getAllCombi()+decalage-1 << "\n";
	for(int i=decalage; i<patterns[indicePattern]->getAllCombi()+decalage ; i++)
	{
		freq[dest][i]=freq[src][i];
	}

//	cerr << "*****Avant****** " << "\n";
//	cerr << "freq[" << index << "][" << premier << "] = " << freq[dest][premier] <<"\n";
	freq[dest][premier] = freq[dest][premier] - 1;
//	cerr << "*****Après****** " << "\n";
//	cerr << "freq[" << index << "][" << premier << "] = " << freq[dest][premier] <<"\n";
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

void FreqKmer::compteFenetre2(int *seq,int seq_taille,int debut,int col,Pattern *p)
{
	int position;
	int i;
	int lk, lp;
	int cpt = 0;
	int flag_premier = 1;
	int indice_kmer = 0;
	int pattern_taille = p->getTaillePattern();
	//~ int flag_inconnu = 0;

	/* On travaille sur une copie, car en param const char */



	/* Parcours de la sequence */
	for(position = debut-1; position < (signed)(seq_taille - pattern_taille); position++)
	{
		indice_kmer = 0;

		/* Pour chaque kmer */
		for(i = position + pattern_taille-1 , lk = 0 ,lp = 0; i >= position; i--, lp++)
		{

			/* si nucléotide inconnu */
			if(seq[i] > 4 || seq[i] < 0)
			{
				/* on va au prochain kmer */
				position = i;

				/* goto bad_nuc : on ne compte pas */
				goto bad_nuc;
			}

			if(p->extraire(pattern_taille-(lp+1)))
			{
				indice_kmer += seq[i] * pow(4,lk);
				lk++;
			}
		}
		cpt++;
		if(flag_premier==1)
		{
			premier= indice_kmer;
			flag_premier = 0;
			cerr << "indice premier toto = " << premier << "\n";
		}
		else
		{
			dernier = indice_kmer;
		}

		freq[index][col+indice_kmer]++;
		bad_nuc: ;

	}
}

void FreqKmer::compteFenetre(int *seq,int seq_taille,int debut ,int indicePattern)
{
	for(int i = debut; i <= debut+seq_taille-patterns[indicePattern]->getTaillePattern() ; i++)
	{
		int col = getCol(indicePattern,seq,i);
		if(i==0)
		{
			premier = col;
		}
		else
		{
			dernier = col;
		}
		freq[index][col]+=1;
	}
}
void FreqKmer::add_one(int *seq,int i,int seq_taille,int indicePattern)
{
	int pattern_taille = patterns[indicePattern]->getTaillePattern();
	int alpha =  i + seq_taille - (pattern_taille + 1);
	int beta = i + seq_taille - 1 ;
	int k;
	int l=0;
	int kmer_taille = patterns[indicePattern]->getTailleKmer();

	int decalage = 0;

	for(int j=0 ; j < indicePattern ; j++)
	{
		decalage += patterns[j]->getAllCombi();
	}

	// printf("seq[%d] = %d\n",beta,seq[beta]);
	// printf("dernier before = %d\n",dernier);
	if (patterns[indicePattern]->isContinue())
	{
		dernier = (dernier-decalage - (seq[alpha]*pow(4,kmer_taille-1)))*4 + seq[beta];
		dernier += decalage;
		//cerr << "dernier = " << dernier << "\n";
	}
	else
	{
		dernier = getCol(indicePattern,seq,i+seq_taille-pattern_taille);
	}

	freq[index][dernier] = freq[index][dernier] + 1;

	premier = 0;
	for (k = 0,l=0; k < pattern_taille; k++)
	{
		if(patterns[indicePattern]->extraire(k))
		{
			premier += seq[i+k] * pow(4,(kmer_taille-l-1));
			l++;
		}
	}
	premier += decalage;
}

void FreqKmer::count(int *seq,int seq_taille,int indicePattern)
{
	int i;
	int taille_sous_sequence = 0;

	premier=0;
	dernier=0;
	if (tailleFenetre==-1 || tailleFenetre>seq_taille)
	{
		taille_sous_sequence = seq_taille;
	}
	else
	{
		taille_sous_sequence = tailleFenetre;
	}

	for(i=0;i<=seq_taille - taille_sous_sequence;i++)
	{
//		cout << "taille_seq = " << seq_taille << "\n";
//		cout << "i = " << i << "\n";
		if (i>0)
		{

			copieLigneFreq(index-1,index,indicePattern);
			add_one(seq,i,taille_sous_sequence,indicePattern);

		}
		else
		{

			//compteFenetre2(seq,taille_sous_sequence,i,col,p);
			compteFenetre(seq,taille_sous_sequence,i,indicePattern);

		}
		index++;
	}
}

void FreqKmer::fillFreq()
{
	initFreq();

	/* Pour chaque fichier fasta */
	for(int k=0;k<nPattern;k++)
	{
		index=0;

		for(int i=0;i<nbFichierFasta;i++)
		{
			/* Pour chaque Data du fichier */
			for(int j=0;j<data[i]->getNtaxa();j++)
			{

				count(data[i]->getDataObject()[j],data[i]->getLengthSeq(j),k);
			}
		}

	}
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

int FreqKmer::getNbLineData(int i)
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

int FreqKmer::getStartLineData(int i)
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

int FreqKmer::getEndLineData(int i)
{
	return indexLineData[i];
}

int FreqKmer::getStartLineDataSeq(int i,int j)
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

int FreqKmer::getEndLineDataSeq(int i,int j)
{
	return indexLineDataSeq[i][j];
}

int FreqKmer::getNbLineDataSeq(int i,int j)
{
	return getEndLineDataSeq(i,j)-getStartLineDataSeq(i,j)+1;
}
