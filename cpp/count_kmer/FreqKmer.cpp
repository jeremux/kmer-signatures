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

using namespace std;

FreqKmer::FreqKmer() {
	/* On veut les accessions donc Yes */
	data=NULL;
	freq=NULL;
	patterns=NULL;
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
	nCol=0;
	nData=0;
	nLigne=0;
	nPattern=0;
	nbFichierFasta=0;
	tailleFenetre=tailleF;
	premier=0;
	dernier=0;
	index=0;

}

FreqKmer::~FreqKmer() {
	delete data;
	for (int var = 0; var < nPattern ; var++)
	{
		delete patterns[var];
	}
	delete[] patterns;
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

	for(int i=0 ; i<nPattern ; i++)
	{

		nCol += patterns[i]->getAllCombi();

	}

}

void FreqKmer::initFromList(string fichier)
{
	string ligne;
	int cpt=0;
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
	ifstream file2(fichier.c_str());
	getline(file2,ligne);
	/* J'init à partir de chaque ligne du fichier */

	while (file2)
	{
		if(tailleLigne!=0)
		{
			data[cpt] = new Data(Yes);
			data[cpt]->initFrom(ligne,Fasta);
			/* Ajouter des attributs ??*/
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
				}
			}
			else
			{
				nLigne +=  data[cpt]->getNtaxa();
			}
			cpt++;
		}
		getline(file2,ligne);
	}
}

void FreqKmer::initFromFasta(string fichier)
{
	nbFichierFasta=1;
	data=new Data*[1];
	data[0] = new Data(Yes);
	data[0]->initFrom(fichier,Fasta);
	nData = data[0]->getNtaxa();
	int tailleSeq=0;
	if (tailleFenetre>0)
	{
		for (int var = 0; var < nData; var++)
		{
			tailleSeq = data[0]->getPrimarySequence(var).length();
			if (tailleSeq<tailleFenetre)
			{
				nLigne += 1;
			}
			else
			{
				nLigne += tailleSeq-tailleFenetre+1;
			}
		}
	}
	else
	{
		nLigne = nData;
	}
}

/**
 * Permet de copier une ligne
 * et mettant a jour les frequences
 * dans le cas d'une graine continue
 */
void FreqKmer::copieLigneFreq(int src, int dest,int col,Pattern *p)
{
	for(int i=col; i<p->getAllCombi()+col ; i++)
	{
		freq[dest][i]=freq[src][i];
	}
	freq[dest][premier+col]--;
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

void FreqKmer::compteFenetre(int *seq,int seq_taille,int debut,int col,Pattern *p)
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
	for(position = debut-1; position < (signed)(seq_taille - pattern_taille ); position++)
	{
		indice_kmer = 0;

		/* Pour chaque kmer */
		for(i = position + pattern_taille , lk = 0 ,lp = 0; i >= position; i--, lp++)
		{

			/* si nucléotide inconnu */
			if(seq[i] > 4 || seq[i] < 0)
			{
				/* on va au prochain kmer */
				position = i;

				/* goto bad_nuc : on ne compte pas */
				goto bad_nuc;
			}


			/* calcul de l'indice */

			if(p->extraire(pattern_taille-(lp+1)))
			{
				// printf("tab_bool[%d] = %d\n", pattern_taille-(lp+1),tab_bool[pattern_taille-(lp+1)]);
				indice_kmer += seq[i] * pow(4,lk);
				lk++;

			}

		}

		/* comptage du kmer */
		cpt++;
		// printf("kmer %d et indice = %d\n",cpt,indice_kmer);
		if(flag_premier)
		{
			premier= indice_kmer;
			flag_premier = 0;
//			cout << "indice premier = " << premier << "\n";
		}

		freq[index][col+indice_kmer]++;
		dernier = indice_kmer;

		bad_nuc: ;

	}

//	cout << "indice dernier = " << dernier << "\n";

}

void FreqKmer::add_one(int *seq,int i,int seq_taille,Pattern *p,int col)
{
	int pattern_taille = p->getTaillePattern();
	int alpha =  i + seq_taille - (pattern_taille + 1);
	int beta = i + seq_taille - 1 ;
	int k;
	int l=0;
	int kmer_taille = p->getTailleKmer();


	// printf("seq[%d] = %d\n",beta,seq[beta]);
	// printf("dernier before = %d\n",dernier);
	if (p->isContinue())
	{
		dernier = (dernier - (seq[alpha]*pow(4,kmer_taille-1)))*4 + seq[beta];
	}
	else
	{
		dernier = 0;
		for (k = i+seq_taille-pattern_taille,l=0; k < i+seq_taille; k++)
		{

			if(p->extraire(k-(i+seq_taille-pattern_taille)))
			{
				dernier += seq[k] * pow(4,(kmer_taille-l-1));
				l++;
				// printf("seq[k] = %d\n",seq[k]);
			}
		}
		// printf("dernier = %d\n",dernier);
	}

	// printf("dernier = %d\n",dernier);
	freq[index][dernier+col] = freq[index][dernier+col] + 1;

	premier = 0;
	for (k = 0,l=0; k < pattern_taille; k++)
	{

		if(p->extraire(k))
		{
			premier += seq[i+k] * pow(4,(kmer_taille-l-1));
			l++;
		}
	}
}

void FreqKmer::count(int *seq,int seq_taille,Pattern *p,int col)
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

			copieLigneFreq(index-1,index,col,p);
			add_one(seq,i,taille_sous_sequence,p,col);
			index++;
		}
		else
		{

			compteFenetre(seq,taille_sous_sequence,i,col,p);
			index++;

		}
	}
}

int FreqKmer::tailleSeq(int *seq,int n)
{
	int res=0;
	for(int i=0;i<n;i++)
	{
		if(seq[i]<4 && seq[i]>=0)
			res++;
	}
	return res;
}

void FreqKmer::fillFreq()
{
	initFreq();
	int col=0;
	/* Pour chaque fichier fasta */
	for(int k=0;k<nPattern;k++)
	{
		index=0;

		for(int i=0;i<nbFichierFasta;i++)
		{
			/* Pour chaque Data du fichier */
			for(int j=0;j<data[i]->getNtaxa();j++)
			{

				count(data[i]->getDataObject()[j],tailleSeq(data[i]->getDataObject()[j],data[i]->getNsite()),patterns[k],col);
			}
		}
		col+=patterns[k]->getAllCombi();
	}
}

void FreqKmer::imprimeCSV()
{
	for(int j=0;j<nCol;j++)
	{
		cout << j << ";";
	}
	cout << "\n";
	for(int i=0;i<nLigne;i++)
	{
		for(int j=0;j<nCol;j++)
		{
			cout << freq[i][j] << ";";
		}
		cout << "\n";
	}
}
