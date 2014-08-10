/*
 * test.cpp
 *
 *  Created on: 17 juil. 2014
 *      Author: jeremy
 */

#include "test.h"
#include "classPattern.h"
#include "classData.h"
#include "FreqKmer.h"

void printResult(bool res,int i)
{
	if (res)
	{
		cerr << "Test" <<i <<" OK !\n";
	}
	else
	{
		cerr << "Test"<<i <<" fail...\n";
	}
}

bool doTest20()
{
	cerr << "\n**test writeCrossVal**\n";
	bool res=true;

	string filename = "Debug/tests/test4/liste.txt";
	string patternPath = "pattern.txt";
	string path_root = "../../create_db/Eukaryota__2759/";
	string key = "genomes";
	cout << "New Freq\n";
	FreqKmer *f = new FreqKmer(-1,patternPath,true,path_root,key);
	cout << "FillFreq\n";
	//FreqKmer *h = f->sampleMe(125);
	//h->fillFreq();
		f->fillFreq();
	f->normalize();
	cout << "writeCrossVal\n";
	f->writeCrossVal(25);

	delete f;
	//delete h;
	printResult(res,20);
	return res;

}

bool doTest19()
{
	cerr << "\n**test init from conf**\n";

	string filename = "Debug/tests/test5/liste.txt";
	string pattern ="Debug/tests/test5/pattern.txt";
	FreqKmer *f = new FreqKmer(5,2,true,filename,pattern,false,"genomes");
	FreqKmer *g=NULL;

	f->fillFreq();

	f->writeConfFeq("Debug/tests/test5/conf.txt");

	g = g->initFromConf("Debug/tests/test5/conf.txt");

	bool res = true;
	res = g->equal(f);
	delete f;
	delete g;
	printResult(res,19);
	return res;
}
bool doTest18()
{
	cerr << "\n**Somme d'une ligne**\n";
	bool res = true;

	string filename = "Debug/tests/test1/liste.txt";
	string pattern = "Debug/tests/test1/pattern.txt";

	FreqKmer *f = new FreqKmer(4,true,filename,pattern,false,"");
	f->fillFreq();
	res = (f->getSum(0,0)==2
		&& f->getSum(0,1)==2
		&& f->getSum(0,2)==2
		&& f->getSum(0,3)==2
		&& f->getSum(0,4)==2
		&& f->getSum(0,5)==2
		&& f->getSum(0,6)==2

		&& f->getSum(1,0)==4
		&& f->getSum(1,1)==4
		&& f->getSum(1,2)==4
		&& f->getSum(1,3)==4
		&& f->getSum(1,4)==4
		&& f->getSum(1,5)==4
		&& f->getSum(1,6)==4
			);

	f->normalize();

	delete f;
	printResult(res,18);
	return res;
}
bool doTest17()
{
	cerr << "\n**Test sample**\n";
	string filename = "Debug/tests/test4/liste.txt";
	FreqKmer *f = new FreqKmer(4,"Debug/tests/test4/pattern.txt",false,"Debug/tests/taxon__alpha","genomes");
	bool res = true;


	FreqKmer *g;
	g = f->sampleMe(1);
	res = res && g->getNbAllTrue()==4;
	delete g;

	g = f->sampleMe(2);
	res = res && g->getNbAllTrue()==8;
	delete g;

	g = f->sampleMe(3);
	res = res && g->getNbAllTrue()==12;
	delete g;

	g = f->sampleMe(4);
	res = res && g->getNbAllTrue()==14;
	delete g;

	g = f->sampleMe(5);
	res = res && g->getNbAllTrue()==16;
	delete g;

	g = f->sampleMe(6);
	res = res && g->getNbAllTrue()==18;
	delete g;

	g = f->sampleMe(7);
	res = res && g->getNbAllTrue()==20;
	delete g;

	g = f->sampleMe(8);
	res = res && g->getNbAllTrue()==22;
	delete g;

	g = f->sampleMe(9);
	res = res && g->getNbAllTrue()==24;
	delete g;

	g = f->sampleMe(10);
	res = res && g->getNbAllTrue()==26;
	delete g;

	g = f->sampleMe(11);
	res = res && g->getNbAllTrue()==27;
	delete g;

	f->fillFreq();
	g = f->sampleMe(18);
	res = res && g->getNbAllTrue()==32;


	g->fillFreq();

	delete g;



	printResult(res,17);
	delete f;


	return res;
}

bool doTest16()
{
	cerr << "\n**Test nombre de sequences a un niveau taxo**\n";
	string filename = "Debug/tests/test4/liste.txt";
	FreqKmer *f = new FreqKmer(4,"Debug/tests/test4/pattern.txt",false,"Debug/tests/taxon__alpha","genomes");
	bool res = true;


	res = (f->getNSeqInTaxa(0) == 3
			&& f->getNSeqInTaxa(1) == 16
			&& f->getNSeqInTaxa(2) == 10
			&& f->getNSeqInTaxa(3) == 3
	);

	printResult(res,16);


	delete f;


	return res;
}

bool doTest15()
{
	cerr << "\n**Test index des chemins fasta dans la list**\n";
	string filename = "Debug/tests/test4/liste.txt";
	FreqKmer *f = new FreqKmer(4,"Debug/tests/test4/pattern.txt",false,"Debug/tests/taxon__alpha","genomes");
	bool res = true;


	res = (f->obtainStartLineTaxaInFastaList(0) == 0
			&& f->obtainStartLineTaxaInFastaList(1) == 1
			&& f->obtainStartLineTaxaInFastaList(2) == 5
			&& f->obtainStartLineTaxaInFastaList(3) == 6

			&& f->obtainEndLineTaxaInFastaList(0) == 0
			&& f->obtainEndLineTaxaInFastaList(1) == 4
			&& f->obtainEndLineTaxaInFastaList(2) == 5
			&& f->obtainEndLineTaxaInFastaList(3) == 7

			&& f->obtainNbLineTaxaInFastaList(0) == 1
			&& f->obtainNbLineTaxaInFastaList(1) == 4
			&& f->obtainNbLineTaxaInFastaList(2) == 1
			&& f->obtainNbLineTaxaInFastaList(3) == 2
	);
	printResult(res,15);

	delete f;


	return res;
}

bool doTest14()
{
	cerr << "\n**test list sous taxon**\n";
	string filename = "Debug/tests/test4/liste.txt";
	FreqKmer *f = new FreqKmer(4,"Debug/tests/test4/pattern.txt",false,"Debug/tests/taxon__alpha","genomes");
	bool res = true;

	res = (f->getPathChildTaxa(2)=="Debug/tests/taxon__alpha/taxon__B"
			&& f->getPathChildTaxa(1)=="Debug/tests/taxon__alpha/taxon__A"
					&& f->getPathChildTaxa(0)=="Debug/tests/taxon__alpha/others"
							&& f->getPathChildTaxa(3)=="Debug/tests/taxon__alpha/taxon__C"
									&& f->getIdTaxa(2)=="B"
											&& f->getIdTaxa(1)=="A"
													&& f->getIdTaxa(0)=="others"
															&& f->getIdTaxa(3)=="C"
																	&& f->getNbChild()==4
	);






	printResult(res,14);

	delete f;

	return res;




}

bool doTest13()
{
	cerr << "\n**test table mask**\n";


	bool res = true;
	string filename = "Debug/tests/test4/liste.txt";
	string pattern = "Debug/tests/test4/pattern.txt";
	string key = "genomes";
	FreqKmer *f = new FreqKmer(4,true,filename,pattern,false,key);

	for(int i=0 ; i < f->getNbFichierFasta(); i++)
	{
		for (int j = 0; j < f->getData()[i]->getNtaxa(); j++)
		{
			res = res && f->takeDataSeq(i,j);
		}
	}

	//f->fillFreq();

	printResult(res,13);
	delete f;

	return res;

}
bool doTest12()
{
	cerr << "\n**test map get DataSeq from line**\n";


	bool res = true;
	string filename = "Debug/tests/test4/liste.txt";
	string pattern = "Debug/tests/test4/pattern.txt";
	string key = "genomes";

	FreqKmer *f = new FreqKmer(4,true,filename,pattern,false,key);


	int idData0,idData1,idData2,idData3,idData4,idData5,idData6,
	idSeq0,idSeq1,idSeq2,idSeq3,idSeq4,idSeq5,idSeq6;

	f->obtainDataSeqFromLine(0,idData0,idSeq0);
	f->obtainDataSeqFromLine(3,idData1,idSeq1);
	f->obtainDataSeqFromLine(5,idData2,idSeq2);
	f->obtainDataSeqFromLine(6,idData3,idSeq3);
	f->obtainDataSeqFromLine(8,idData4,idSeq4);
	f->obtainDataSeqFromLine(11,idData5,idSeq5);
	f->obtainDataSeqFromLine(14,idData6,idSeq6);









	res = (idData0 == 0)
		&& (idData1 == 0)
		&& (idData2 == 1)
		&& (idData3 == 1)
		&& (idData4 == 1)
		&& (idData5 == 1)
		&& (idData6 == 1)
		&& (idSeq0 == 0)
		&& (idSeq1 == 1)
		&& (idSeq2 == 0)
		&& (idSeq3 == 0)
		&& (idSeq4 == 1)
		&& (idSeq5 == 2)
		&& (idSeq6 == 2);

	//f->fillFreq();


	printResult(res,12);


	delete f;
	return res;
}

bool doTest11()
{
	cerr << "\n**Données réelles validées par perl/C**\n";


	string filename = "Debug/tests/haemo/haemo.fasta";
	string pattern = "Debug/tests/haemo/pattern5.txt";
	FreqKmer *f = new FreqKmer(-1,1,false,filename,pattern,false,"");

	int val[] = {11,10,7,18,12,5,5,7,10,3,7,6,29,4,14,17,8,7,9,14,7,2,3,6,4,2,3,5,4,6,8,13,14,4,9,7,2,3,1,5,5,1,0,9,7,1,3,9,26,14,8,31,4,8,1,4,10,8,6,19,24,8,13,22,8,5,7,10,9,3,1,4,7,5,7,4,16,5,7,19,9,4,4,11,
			4,2,1,1,1,3,1,9,7,4,3,10,5,1,4,3,1,2,0,4,0,0,0,3,6,1,1,4,10,7,4,12,8,1,2,6,1,3,6,6,3,8,5,13,9,7,7,11,2,2,4,2,6,1,2,3,13,4,10,6,5,3,6,6,7,1,3,6,1,0,1,2,7,5,4,4,8,4,5,10,2,1,2,1,2,1,0,2,5,1,2,8,5,2,7,19,2,1,2,1,3,2,4,4,11,7,5,9,28,8,3,35,16,11,3,
			15,10,10,6,11,34,11,16,46,8,3,6,13,7,3,1,10,2,4,1,1,1,5,6,7,6,4,7,15,13,6,1,14,20,2,8,7,21,8,5,20,27,
			23,17,38,15,8,2,22,8,8,8,7,31,19,18,25,9,3,8,12,7,1,1,10,13,5,3,3,11,4,9,13,4,4,0,8,3,0,3,4,2,1,0,2,5,
			2,1,4,8,3,0,9,5,6,0,5,8,4,2,3,9,2,3,3,7,9,8,20,6,4,3,4,4,10,8,9,16,14,3,19,8,4,8,8,3,3,0,4,5,3,5,5,6,3,
			9,13,5,0,5,3,2,0,0,1,1,0,0,0,0,0,2,3,0,2,0,0,0,2,0,2,1,0,2,1,5,3,0,4,2,4,1,12,3,2,1,0,7,0,4,4,3,2,4,
			12,0,4,0,4,3,0,0,3,0,3,0,2,2,1,1,2,1,2,2,1,1,1,2,2,0,0,0,1,1,1,3,4,2,0,1,0,0,0,0,1,0,2,1,0,3,1,1,3,5,
			3,2,5,2,1,1,2,1,0,0,1,6,1,1,7,12,5,1,6,8,7,1,5,1,5,2,2,15,1,9,15,7,3,2,12,4,1,0,1,2,1,1,3,7,5,0,4,6,2,
			2,6,1,0,1,2,3,2,3,8,7,6,3,6,5,5,2,11,7,7,2,6,2,0,2,6,16,6,8,15,8,6,1,7,7,1,2,9,3,1,2,8,14,3,2,9,5,0,1,
			3,1,1,4,2,1,4,0,2,4,2,3,4,2,3,1,5,3,3,0,3,4,0,0,1,5,2,3,8,11,8,6,16,3,3,1,3,4,7,5,11,14,5,3,11,5,4,5,6,
			2,2,0,2,5,3,1,3,7,3,3,6,4,0,3,7,2,0,0,1,0,1,3,2,6,1,2,3,1,1,0,1,0,0,0,0,0,0,0,2,1,0,1,3,4,4,1,7,6,2,1,3,
			2,0,3,5,7,6,0,9,6,4,3,7,2,1,1,5,3,4,0,6,14,1,9,16,2,0,1,2,2,1,0,1,0,0,1,2,3,1,1,1,3,3,1,4,1,2,1,1,3,1,1,1,2,2,1,4,8,2,2,9,3,1,2,3,1,4,4,3,6,7,4,7,9,6,3,12,9,3,0,6,8,0,2,
			8,13,7,9,22,4,1,3,6,4,0,0,1,1,1,2,3,1,3,1,4,2,3,2,2,3,1,0,3,5,0,1,5,6,0,0,6,9,5,3,18,13,5,3,6,6,3,4,3,16,7,3,7,18,11,10,27,12,11,7,7,8,2,3,3,28,4,18,27,13,6,13,23,17,5,4,12,7,
			0,0,3,21,7,4,8,10,0,2,12,10,5,3,7,10,1,3,3,12,1,3,12,30,14,13,40,17,7,2,8,14,8,18,15,51,20,12,40,11,6,5,
			14,1,2,4,3,2,5,4,
			7,16,7,11,14,10,6,7,10,5,1,0,2,0,0,0,1,6,1,9,6,1,2,1,2,5,2,1,3,2,1,1,2,4,2,0,4,9,6,4,11,7,1,3,7,6,0,4,6,10,6,2,11,7,4,4,6,2,5,2,3,2,1,3,7,12,4,6,9,12,1,3,10,4,0,1,3,2,0,0,0,5,5,3,12,7,2,6,26,2,1,0,3,6,1,4,6,11,5,7,10,11,10,7,18,7,2,2,4,3,1,3,4,12,12,6,11,16,17,9,24,21,17,6,14,5,10,7,7,35,14,
			21,39,16,3,7,17,18,4,0,10,1,5,2,3,20,5,9,14,7,3,2,8,9,1,0,6,13,2,5,13,13,0,3,9,24,25,7,41,9,11,4,
			14,4,4,19,9,34,6,7,14
	};

	int ligne = 1;
	int col = 1024;
	double ** freq = new double *[1];
	freq[0] = new double[1024];
	for(int j=0 ; j<col ; j++)
	{
		freq[0][j]=val[j];
	}

	f->fillFreq();
	if (f->getNLine()!=1 || f->getNCol()!=1024)
	{
		cerr << "Test sur nLigne non ok\n";
	}


	//      cout << "Fin fillFreq\n";
	//      cout << "Debut Impression\n";
	bool res=true;
	for(int i=0;i<ligne;i++)
	{
		for(int j=0;j<col;j++)
		{
			//                		cerr << "freq[" << i << "][" << j << "] = " << f->getFreq()[i][j] <<  " || fb[" << i << "][" << j << "] = " << freq[i][j]<< "\n";
			if(f->getFreq()[i][j]!=freq[i][j])
			{
				res=false;

			}
			if(!res)
				break;
		}
		if(!res)
			break;
	}
	printResult(res,11);



	for(int i=0 ; i<ligne ; i++)
	{
		delete[] freq[i];
	}
	delete[] freq;

	delete f;
	return res;
}

bool doTest10()
{
	cerr << "\n**test decalage bis**\n";


	string filename = "Debug/tests/test1/liste2.txt";
	string pattern = "Debug/tests/test1/pattern2.txt";

	FreqKmer *f = new FreqKmer(3,2,true,filename,pattern,false,"");




	int ligne = 2;
	int col = 16;
	double ** freq = new double *[ligne];
	for(int i=0 ; i<ligne ; i++)
	{
		freq[i] = new double[col];
		for(int j=0 ; j<col ; j++)
		{
			freq[i][j]=0;
		}
	}

	freq[0][0] = 1;

	freq[1][3] = 1;


	f->fillFreq();

	if (f->getNLine()!=2 || f->getNCol()!=16)
	{
		cerr << "Test sur nLigne non ok\n";
	}
	//      cout << "Fin fillFreq\n";
	//      cout << "Debut Impression\n";
	bool res=true;
	for(int i=0;i<ligne;i++)
	{
		for(int j=0;j<col;j++)
		{
			// cerr << "freq[" << i << "][" << j << "] = " << f->getFreq()[i][j] <<  " || fb[" << i << "][" << j << "] = " << freq[i][j]<< "\n";
			if(f->getFreq()[i][j]!=freq[i][j])
			{
				res=false;

			}
			if(!res)
				break;
		}
		if(!res)
			break;
	}
	printResult(res,10);



	for(int i=0 ; i<ligne ; i++)
	{
		delete[] freq[i];
	}
	delete[] freq;

	delete f;
	return res;
}

bool doTest9()
{
	cerr << "\n**test decalage 2**\n";

	string filename = "Debug/tests/test5/liste.txt";
	string pattern ="Debug/tests/test5/pattern.txt";
	FreqKmer *o = new FreqKmer(5,2,true,filename,pattern,false,"");



	int ligne = 5;
	int col = 20;
	double ** freq = new double *[ligne];
	for(int i=0 ; i<ligne ; i++)
	{
		freq[i] = new double[col];
		for(int j=0 ; j<col ; j++)
		{
			freq[i][j]=0.0;
		}
	}

	freq[0][0] = 1 ;
	freq[0][3] = 1 ;
	freq[0][5] = 1 ;
	freq[0][16] = 2 ;
	freq[0][17] = 2 ;
	freq[0][19] = 1 ;
	freq[1][3] = 1 ;
	freq[1][6] = 1 ;
	freq[1][12] = 1 ;
	freq[1][16] = 2 ;
	freq[1][17] = 1 ;
	freq[1][18] = 1 ;
	freq[1][19] = 1 ;
	freq[2][3] = 1 ;
	freq[2][12] = 1 ;
	freq[2][13] = 1 ;
	freq[2][16] = 2 ;
	freq[2][17] = 1 ;
	freq[2][19] = 2 ;
	freq[3][1] = 1 ;
	freq[3][6] = 1 ;
	freq[3][12] = 1 ;
	freq[3][16] = 1 ;
	freq[3][17] = 2 ;
	freq[3][18] = 1 ;
	freq[3][19] = 1 ;
	freq[4][1] = 1 ;
	freq[4][7] = 1 ;
	freq[4][8] = 1 ;
	freq[4][16] = 2 ;
	freq[4][17] = 1 ;
	freq[4][18] = 1 ;
	freq[4][19] = 1 ;


	o->fillFreq();

	if (o->getNLine()!=5 || o->getNCol()!=20)
	{
		cerr << "Test sur nLigne non ok\n";
	}
	//      cout << "Fin fillFreq\n";
	//      cout << "Debut Impression\n";
	bool res=true;
	for(int i=0;i<ligne;i++)
	{
		for(int j=0;j<col;j++)
		{
			if(o->getFreq()[i][j]!=freq[i][j])
			{
				res=false;
				break;
			}
		}
	}
	printResult(res,9);


	for(int i=0 ; i<ligne ; i++)
	{
		delete[] freq[i];
	}
	delete[] freq;

	delete o;
	return res;
}

bool doTest8()
{
	cerr << "\n**test decalage**\n";

	string filename = "Debug/tests/test1/liste.txt";
	string pattern ="Debug/tests/test1/pattern.txt";
	FreqKmer *f = new FreqKmer(3,2,true,filename,pattern,false,"");



	unsigned int ligne = 5;
	unsigned int col = 20;

	double **freq = new double *[ligne];
	for(unsigned int i=0 ; i<ligne ; i++)
	{
		freq[i] = new double[col];
		for(unsigned int j=0 ; j<col ; j++)
		{
			freq[i][j]=0;
		}
	}

	freq[0][0] = 1;

	freq[1][3] = 1;

	freq[2][3] = 1;

	freq[3][12] = 1;

	freq[4][1] = 1;


	freq[0][16] = 2;
	freq[0][17] = 1;

	freq[1][16] = 1;
	freq[1][17] = 1;
	freq[1][19] = 1;

	freq[2][16] = 1;
	freq[2][19] = 2;

	freq[3][16] = 1;
	freq[3][17] = 1;
	freq[3][19] = 1;

	freq[4][16] = 1;
	freq[4][17] = 1;
	freq[4][18] = 1;


	f->fillFreq();

	if (f->getNLine()!=5 || f->getNCol()!=20)
	{
		cerr << "Test sur nLigne non ok\n";
	}
	//      cout << "Fin fillFreq\n";
	//      cout << "Debut Impression\n";
	bool res=true;
	for(unsigned int i=0;i<ligne;i++)
	{
		for(unsigned int j=0;j<col;j++)
		{
			//                		cerr << "freq[" << i << "][" << j << "] = " << f->getFreq()[i][j] <<  " || fb[" << i << "][" << j << "] = " << freq[i][j]<< "\n";
			if(f->getFreq()[i][j]!=freq[i][j])
			{
				res=false;

			}
			if(!res)
				break;
		}
		if(!res)
			break;
	}
	printResult(res,8);

	for(unsigned int i=0 ; i<ligne ; i++)
	{
		delete[] freq[i];
	}


	delete[] freq;
	delete f;
	return res;
}

bool doTest7()
{
	cerr << "\n**Test de la map deux**\n";

	bool res = true;
	string filename = "Debug/tests/test4/liste.txt";
	string pattern = "Debug/tests/test4/pattern.txt";
	FreqKmer *f = new FreqKmer(4,true,filename,pattern,false,"");

	res = (
			f->obtainStartLineDataSeq(0,0)==0
			&& f->obtainEndLineDataSeq(0,0)==2
			&& f->obtainNbLineDataSeq(0,0)==3

			&& f->obtainStartLineDataSeq(0,1)==3
			&& f->obtainEndLineDataSeq(0,1)==3
			&& f->obtainNbLineDataSeq(0,1)==1

			&& f->obtainStartLineDataSeq(1,0)==4
			&& f->obtainEndLineDataSeq(1,0)==7
			&& f->obtainNbLineDataSeq(1,0)==4

			&& f->obtainStartLineDataSeq(1,1)==8
			&& f->obtainEndLineDataSeq(1,1)==9
			&& f->obtainNbLineDataSeq(1,1)==2

			&& f->obtainStartLineDataSeq(1,2)==10
			&& f->obtainEndLineDataSeq(1,2)==14
			&& f->obtainNbLineDataSeq(1,2)==5
	);



	printResult(res,7);


	delete f;
	return res;
}

bool doTest6()
{
	cerr << "\n**Test de la map une**\n";

	bool res = true;
	string filename = "Debug/tests/test1/liste.txt";
	string pattern = "Debug/tests/test1/pattern.txt";
	FreqKmer *f = new FreqKmer(4,true,filename,pattern,false,"");

	res = !(f->obtainNbLineData(0)!=3 || f->obtainNbLineData(1)!=4 || f->obtainStartLineData(0)!=0 || f->obtainEndLineData(0)!=2 || f->obtainStartLineData(1)!=3 || f->obtainEndLineData(1)!=6);

	printResult(res,6);


	delete f;
	return res;
}
bool doTest5()
{
	cerr << "\n**Test espace des kmers**\n";
	bool res = true;
	string filename = "Debug/tests/test1/liste.txt";
	string pattern = "Debug/tests/test1/pattern.txt";

	FreqKmer *f = new FreqKmer(4,true,filename,pattern,false,"");


	if (f->obtainStartColKmer(0)!=0 || f->obtainStartColKmer(1)!=16 || f->obtainEndColKmer(0)!=15 || f->obtainEndColKmer(1)!=19)
	{
		res=false;
	}

	printResult(res,5);


	delete f;
	return res;
}
bool doTest4()
{
	cerr << "\n**Test taille des sequences**\n";

	bool res = true;
	string filename = "Debug/tests/test1/liste.txt";
	string pattern = "Debug/tests/test1/pattern.txt";
	FreqKmer *f = new FreqKmer(4,true,filename,pattern,false,"");

	if (f->getData()[0]->getLengthSeq(0)!=6 || f->getData()[1]->getLengthSeq(0)!=7)
	{
		res=false;
	}

	printResult(res,4);


	delete f;
	return res;
}
bool doTest3()
{
	cerr << "\n**test multi-pattern**\n";


	string filename = "Debug/tests/test1/liste.txt";
	string pattern = "Debug/tests/test1/pattern.txt";

	FreqKmer *f = new FreqKmer(4,true,filename,pattern,false,"");

	int ligne = 7;
	int col = 20;
	double ** freq = new double *[ligne];
	for(int i=0 ; i<ligne ; i++)
	{
		freq[i] = new double[col];
		for(int j=0 ; j<col ; j++)
		{
			freq[i][j]=0.0;
		}
	}
	freq[0][0] = 1;
	freq[0][5] = 1;
	freq[1][3] = 1;
	freq[1][5] = 1;
	freq[2][3] = 1;
	freq[2][6] = 1;
	freq[3][3] = 1;
	freq[3][13] = 1;
	freq[4][12] = 1;
	freq[4][13] = 1;
	freq[5][6] = 1;
	freq[5][12] = 1;
	freq[6][6] = 1;
	freq[6][1] = 1;
	freq[0][16] = 2;
	freq[0][17] = 2;
	freq[1][16] = 1;
	freq[1][17] = 2;
	freq[1][19] = 1;
	freq[2][16] = 1;
	freq[2][17] = 1;
	freq[2][18] = 1;
	freq[2][19] = 1;
	freq[3][16] = 1;
	freq[3][17] = 1;
	freq[3][19] = 2;
	freq[4][16] = 1;
	freq[4][17] = 1;
	freq[4][19] = 2;
	freq[5][16] = 1;
	freq[5][17] = 1;
	freq[5][18] = 1;
	freq[5][19] = 1;
	freq[6][16] = 1;
	freq[6][17] = 2;
	freq[6][18] = 1;

	f->fillFreq();
	//      cout << "Fin fillFreq\n";
	//      cout << "Debut Impression\n";
	bool res=true;
	for(int i=0;i<ligne;i++)
	{
		for(int j=0;j<col;j++)
		{
			if(f->getFreq()[i][j]!=freq[i][j])
			{
				res=false;
				break;
			}
		}
	}
	if (res)
	{
		cerr << "Test3 OK !\n";
	}
	else
	{
		cerr << "Test3 fail...\n";
	}

	for(int i=0 ; i<ligne ; i++)
	{
		delete[] freq[i];
	}


	delete[] freq;

	delete f;
	return res;
}

bool doTest2()
{
	cerr << "\n**test pattern simple**\n";

	string filename = "Debug/tests/test1/liste.txt";
	string pattern = "Debug/tests/test2/pattern2.txt";

	FreqKmer *g = new FreqKmer(4,true,filename,pattern,false,"");


	int ligne = 7;
	int col = 16;
	double ** freq = new double *[ligne];
	for(int i=0 ; i<ligne ; i++)
	{
		freq[i] = new double[col];
		for(int j=0 ; j<col ; j++)
		{
			freq[i][j]=0.0;
		}
	}
	freq[0][0] = 1;
	freq[0][5] = 1;
	freq[1][3] = 1;
	freq[1][5] = 1;
	freq[2][3] = 1;
	freq[2][6] = 1;
	freq[3][3] = 1;
	freq[3][13] = 1;
	freq[4][12] = 1;
	freq[4][13] = 1;
	freq[5][6] = 1;
	freq[5][12] = 1;
	freq[6][6] = 1;
	freq[6][1] = 1;

	g->fillFreq();


	//      cout << "Debut Impression\n";
	bool res=true;
	for(int i=0;i<ligne;i++)
	{
		for(int j=0;j<col;j++)
		{
			if(g->getFreq()[i][j]!=freq[i][j])
			{
				res=false;
				break;
			}
		}
	}
	if (res)
	{
		cerr << "Test2 OK !\n";
	}
	else
	{
		cerr << "Test2 fail...\n";
	}
	for(int i=0 ; i<ligne ; i++)
	{
		delete[] freq[i];
	}

	delete[] freq;
	delete g;
	return res;
}

bool doTest1()
{
	cerr << "\n**test dimension tableau**\n";

	string filename = "Debug/tests/test3/test.fasta";

	string pattern = "Debug/tests/test3/3patterns.txt";

	FreqKmer *f = new FreqKmer(6,false,filename,pattern,false,"");
	bool res = true;

	if(f->getNCol()!=84 || f->getNLine()!=107)
	{
		cerr << "Test1 fail...\n";
		res = false;
	}
	else
	{
		cerr << "Test1 OK !\n";
	}
	f->initFreq();

	delete f;
	return res;
}

bool doTest0()
{
	cerr << "\n**test new Data**\n";
	Switch s;
	bool res = false;
	Data *d1 = new Data();


	s = d1->initFrom("Debug/tests/test1/file3.fasta",Fasta);
	if(s==Yes)
		res=true;

	printResult(res,0);

	delete d1;
	return res;


}
bool callTest(int i)
{
	bool res = false;
	printSwitch(getSwitch(0));
	switch (i) {
	case 0:
		res = doTest0();
		break;

	case 1:
		res = doTest1();
		break;
	case 2:
		res = doTest2();
		break;
	case 3:
		res = doTest3();
		break;
	case 4:
		res = doTest4();
		break;
	case 5:
		res = doTest5();
		break;
	case 6:
		res = doTest6();
		break;
	case 7:
		res = doTest7();
		break;
	case 8:
		res = doTest8();
		break;
	case 9:
		res = doTest9();
		break;
	case 10:
		res = doTest10();
		break;
	case 11:
		res = doTest11();
		break;
	case 12:
		res = doTest12();
		break;
	case 13:
		res = doTest13();
		break;
	case 14:
		res = doTest14();
		break;
	case 15:
		res = doTest15();
		break;
	case 16:
		res = doTest16();
		break;
	case 17:
		res = doTest17();
		break;
	case 18:
		res = doTest18();
		break;
	case 19:
		res = doTest19();
		break;
	case 20:
		res = doTest20();
		break;
	default:
		break;
	}
	return res;
}



bool executeTests(int k)
{
	bool res = true ;
	cerr << "===Début des tests===\n";
	for(int i=0; i <= k ; i++)
	{
		res = res & callTest(i);
	}
	cerr << "\n===Fin des tests===\n";
	return res;

}


