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

/******************************************************************************
 *
 *Constructeur
 ****************************************************************************/


FreqKmer::FreqKmer(int win_size,int s,bool list, string file,string patternFile, bool b,string key)
{

	printSwitch(getSwitch(0));
	key_fasta=key;
	indexTaxaInFasta=NULL;
	freqFilled=false;
	pathRoot="null";
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

	}




}


/******************************************************************************
 *
 *Constructeur
 ****************************************************************************/
FreqKmer::FreqKmer(int win_size,bool list, string file,string patternFile, bool b,string key)
{

	printSwitch(getSwitch(0));
	key_fasta=key;
	freqFilled=false;
	mask=NULL;
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


FreqKmer::FreqKmer(int win_size,string patternFile, bool b,string pathR,string key)
{

	printSwitch(getSwitch(0));
	key_fasta=key;
	pathRoot=pathR;
	freqFilled=false;
	mask=NULL;
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

FreqKmer::FreqKmer(int win_size,int s,string patternFile, bool b,string pathR,string key)
{

	printSwitch(getSwitch(0));
	key_fasta=key;
	pathRoot=pathR;
	freqFilled=false;
	mask=NULL;
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





	if(dataVerbose){
		cerr << "Debut FreqKmer::~FreqKmer()\n ";
		cerr.flush();
	}
}

/******************************************************************************
 *
 *Initialise le fichier de pattern
 ****************************************************************************/
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

/******************************************************************************
 *
 *Constructeur
 ****************************************************************************/
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

	mask[0] = new bool[nSeq];

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
 * @param	debut: où commencer à compter
 * @param	indexPattern: indice du kmer courant
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

			if(freq[start_line]==NULL)
			{
				//				cout << "ALLOC freq[" << start_line << "]\n";
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
 * @param 	seq: la sequence ou compter
 * @param 	seq_length : taille de la sequence
 * @parm	indexPattern : indice du kmer actuel
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

		nTaxa = data[i]->getNtaxa();

		if(noData)
		{
			delete data[i];
			data[i]=NULL;

		}

		for(int j=0;j<nTaxa;j++)
		{

			fillFreq(i,j);
		}



	}


	freqFilled=true;

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
		//		myfile << j << ";";
		//		cout << j << ";";
	}
	//	myfile << endl;
	//	cout << endl;
	int indexData,indexSeq;
	for(int i=0;i<nLine;i++)
	{
		obtainDataSeqFromLine(i,indexData,indexSeq);
		if(mask[indexData][indexSeq])
		{

			for(int j=0;j<nCol;j++)

			{
				myfile << freq[i][j] << ";";
			}
			myfile << endl;
		}
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

bool FreqKmer::takeDataSeq(int i,int j)
{
	return mask[i][j];
}

void FreqKmer::setFalseMask(int i,int j)
{
	mask[i][j] = false;
}
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
	while ((dirp = readdir(dp)) != NULL)
	{
		s1 = string(dirp->d_name);

		if (s1.find("__") != std::string::npos || s1=="others")
		{
			files.push_back(dir+"/"+s1);


		}

	}
	sort(files.begin(),files.end());

	for(unsigned i=0;i< files.size();i++)
	{
		s1 = files[i] ;
		found=s1.find_last_of("/\\");
		if(found==s1.size()-1)
		{
			// on enlève le slash si en fin de path : ex taxonA/taxonB/taxonC/ "
			// ensuite vaut taxonA/taxonB/taxonC
			s1 = s1.substr(found,1);
			found=s1.find_last_of("/\\");
		}
		s1=s1.substr(found+1,s1.size());
		if (s1=="others")
		{
			taxids.push_back(s1);
		}
		else
		{
			x = s1.find( "__" ) + 2 ;
			taxids.push_back(s1.substr(x,s1.length()));
		}

	}
	nbChildTaxa = files.size();
	closedir(dp);
	return 0;
}

string FreqKmer::getIdTaxa(int dir_i)
{
	return idTaxa[dir_i];
}

string FreqKmer::getPathChildTaxa(int dir_i)
{
	return pathChildTaxa[dir_i];
}

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

int FreqKmer::obtainEndLineTaxaInFastaList(int i)
{
	return indexTaxaInFasta[i];
}


int FreqKmer::obtainNbLineTaxaInFastaList(int i)
{
	int res = 0;
	res = (obtainEndLineTaxaInFastaList(i)-obtainStartLineTaxaInFastaList(i))+1;

	//	cerr << "obtainEndLineTaxaInFastaList(" << i << ") = " << obtainEndLineTaxaInFastaList(i) << "\n";
	//	cerr << "obtainStartLineTaxaInFastaList(" << i << ") = " << obtainStartLineTaxaInFastaList(i) << "\n";
	//	cerr << "Pour i = " << i << " je retourne " << res << "\n";
	return res;
}

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
	for(int i=0; i < nbChildTaxa ; i++)
	{
		path = getPathChildTaxa(i);

		path += "/data/fasta/nucleotides/" + key_fasta;

		if((dp  = opendir(path.c_str()) )== NULL)
		{
			cerr << "Error(" << errno << ") opening " << path << endl;
			exit(0);
		}
		while ((dirp = readdir(dp)) != NULL)
		{
			file_name=dirp->d_name;
			path_tmp = path+"/"+file_name;

			if(!(dirp->d_type == DT_DIR))
			{

				x = file_name.find( "." ) + 1 ;
				extension = file_name.substr(x,file_name.length());
				if(extension=="fasta")
				{

					listPathFasta.push_back(path_tmp);

					cpt++;
				}

			}
		}

		//	indexTaxaInFasta[i]=cpt;
		closedir(dp);
	}
	sort(listPathFasta.begin(),listPathFasta.end());

	ancien=listPathFasta[0];

	found=ancien.find_last_of("/\\");
	previous=ancien.substr(0,found);

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
		//		cerr << "Jécris " << listPathFasta[j] << "\n";
		myfile << listPathFasta[j] << "\n";

	}
	myfile.close();
	if (dataVerbose)
	{
		cerr << "Fin FreqKmer::writeListFasta()\n";
	}
}

void FreqKmer::randomTab(vector<int> *result,int tabSize,int sampleSize)
{
	if(dataVerbose)
	{
		cerr << "DEBUT FreqKmer::randomTab(tab,"<<tabSize<<","<<sampleSize<<")\n";
	}
	result->erase(result->begin(),result->end());
	if(sampleSize>tabSize)
	{
		return;
	}

	int *tmp = new int[tabSize];
	int r = -1;
	int sup = tabSize;
	int val_tmp;



	for(int i=0;i<tabSize;i++)
	{
		tmp[i]=i;
	}


	for(int j=0;j<sampleSize;j++)
	{

		r = rand() % sup;
		val_tmp = tmp[sup-1];
		result->push_back(tmp[r]);
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

FreqKmer* FreqKmer::sampleMe(int sampleSize)
{
	FreqKmer *res;

	sampledTaxon.erase(sampledTaxon.begin(),sampledTaxon.end());

	if (dataVerbose)
	{
		cerr << "Debut FreqKmer::sampleMe("<< sampleSize << ")\n";
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
		if(noData)
		{
			if(data[i]!=NULL)
			{
				cerr << "WARNING in FreqKmer::sampleMe, NULL pointer on data[" << i << "] expected\n";
				exit(0);
			}
			data[i] = new Data();
			data[i]->initFrom(pathFasta[i],Fasta);
		}
		nbSequences = data[i]->getNtaxa();
		mask_tmp[i] = new bool[nbSequences];
		for(int j=0;j<nbSequences;j++)
		{
			mask_tmp[i][j]=false;
		}

		if(noData)
		{
			delete data[i];
			data[i]=NULL;
		}
	}

	for(int i=0;i<nbChildTaxa;i++)
	{
		nbSeqTaxa = getNSeqInTaxa(i);
		/* si on doit tirer plus de qu'il y a
		 * de seq alors on tire tout
		 */

		if(sampleSize>=nbSeqTaxa)
		{

			d = obtainStartLineTaxaInFastaList(i);
			f = obtainEndLineTaxaInFastaList(i);

			for(int j=d; j<=f;j++)
			{
				if(noData)
				{
					if(data[j]!=NULL)
					{
						cerr << "WARNING in FreqKmer::sampleMe, NULL pointer on data[" << j << "] expected\n";
						exit(0);
					}
					data[j] = new Data();
					data[j]->initFrom(pathFasta[j],Fasta);
				}
				for(int k=0;k<data[j]->getNtaxa();k++)
				{
					mask_tmp[j][k]=true;
					sampledTaxon.push_back(intPair(j,k));

				}
				if(noData)
				{
					delete data[j];
					data[j]=NULL;
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
		if(noData)
		{
			if(data[i]!=NULL)
			{
				cerr << "WARNING in FreqKmer::sampleMe, NULL pointer on data[" << i << "] expected\n";
				exit(0);
			}
			data[i] = new Data();
			data[i]->initFrom(pathFasta[i],Fasta);
		}
		nbSequences = data[i]->getNtaxa();
		if(noData)
		{
			delete data[i];
			data[i]=NULL;
		}
		for(int j=0;j<nbSequences;j++)
		{
			res->mask[i][j]=mask_tmp[i][j];
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
	return res;

}

int FreqKmer::getNSeqInTaxa(int i)
{
	int nbData = 0;
	nbData = obtainNbLineTaxaInFastaList(i) ;
	//	cerr << "nbData = " << nbData << "\n";
	int res=0;
	int startLineInFasta = 0;
	startLineInFasta = obtainStartLineTaxaInFastaList(i);
	for(int j=0; j<nbData ; j++)
	{
		if(noData)
		{
			if(data[j+startLineInFasta]!=NULL)
			{
				cerr << "WARNING in FreqKmer::sampleMe, NULL pointer on data[" << j+startLineInFasta << "] expected\n";
				exit(0);
			}
			data[j+startLineInFasta] = new Data();
			data[j+startLineInFasta]->initFrom(pathFasta[j+startLineInFasta],Fasta);
		}
		res+=data[j+startLineInFasta]->getNtaxa();
		if(noData)
		{
			delete data[j+startLineInFasta];
			data[j+startLineInFasta]=NULL;
		}
	}

	return res;
}

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

int FreqKmer::getNbAllTrue()
{
	int res = 0;
	for(int i=0;i<nbFastaFile;i++)
	{
		if(noData)
		{
			if(data[i]!=NULL)
			{
				cerr << "WARNING in FreqKmer::sampleMe, NULL pointer on data[" << i << "] expected\n";
				exit(0);
			}
			data[i] = new Data();
			data[i]->initFrom(pathFasta[i],Fasta);
		}
		res+=getNbTrue(mask[i],0,data[i]->getNtaxa()-1);
		if(noData)
		{
			delete data[i];
			data[i]=NULL;
		}
	}
	return res;
}

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
		if(noData)
		{
			if(data[i]!=NULL)
			{
				cerr << "WARNING in FreqKmer::sampleMe, NULL pointer on data[" << i << "] expected\n";
				exit(0);
			}
			data[i] = new Data();
			data[i]->initFrom(pathFasta[i],Fasta);
		}
		nSeq = data[i]->getNtaxa();
		if(noData)
		{
			delete data[i];
			data[i]=NULL;
		}

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
	int nbComb = patterns[indexPattern]->getAllCombi();

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
		if(noData)

		{
			if(data[i]!=NULL)
			{
				cerr << "WARNING in FreqKmer::writeCrossVal, NULL pointer on data[" << i << "] expected\n";
				exit(0);
			}
			data[i] = new Data();
			data[i]->initFrom(pathFasta[i],Fasta);
		}
		seqInData_i= data[i]->getNtaxa();
		if(noData)
		{
			delete data[i];
			data[i]=NULL;
		}
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
		for(int t=0;t<nPattern;t++)
		{
			normalizeLine(t,l);
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
		cerr << "FreqKmer::fillFreq()\n ";
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
		cerr << "Fin FreqKmer::fillFreq()\n ";
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
		if(noData)
		{
			if(data[i]!=NULL)
			{
				cerr << "WARNING in FreqKmer::sampleMe, NULL pointer on data[" << i << "] expected\n";
				exit(0);
			}
			data[i] = new Data();
			data[i]->initFrom(pathFasta[i],Fasta);
		}
		nbSequences = data[i]->getNtaxa();

		for(int j=0;j<nbSequences;j++)
		{
			res->mask[i][j]=false;

		}



		if(noData)
		{
			delete data[i];
			data[i]=NULL;
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
		if(noData)

		{
			if(data[i]!=NULL)
			{
				cerr << "WARNING in FreqKmer::writeCrossVal, NULL pointer on data[" << i << "] expected\n";
				exit(0);
			}
			data[i] = new Data();
			data[i]->initFrom(pathFasta[i],Fasta);
		}
		seqInData_i= data[i]->getNtaxa();
		if(noData)
		{
			delete data[i];
			data[i]=NULL;
		}
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
		if(noData)

		{
			if(data[i]!=NULL)
			{
				cerr << "WARNING in FreqKmer::writeCrossVal, NULL pointer on data[" << i << "] expected\n";
				exit(0);
			}
			data[i] = new Data();
			data[i]->initFrom(pathFasta[i],Fasta);
		}
		seqInData_i= data[i]->getNtaxa();
		if(noData)
		{
			delete data[i];
			data[i]=NULL;
		}
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
		for(int t=0;t<f->nPattern;t++)
		{
			f->normalizeLine(t,l);
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
		toLearn = this;
		toPredict = toLearn;
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
	delete toLearn;

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
			toPredict = toLearn;
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

