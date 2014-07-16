/*
 * count_kmer.cpp
 *
 *  Created on: 10 juil. 2014
 *      Author: jeremy
 */
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <typeinfo>
#include <iomanip>


#include "classPattern.h"
#include "classData.h"
#include "FreqKmer.h"

void doTest3()
{
	cerr << "\n**test multi-pattern**\n";

	FreqKmer *f = new FreqKmer(4);
	string filename = "tests/test1/liste.txt";
	f->initPatterns("tests/test1/pattern.txt");
	f->initFromList(filename);

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
//	cout << "Fin fillFreq\n";
//	cout << "Debut Impression\n";
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
	delete f;
}

void doTest2()
{
	cerr << "\n**test pattern simple**\n";
	FreqKmer *g = new FreqKmer(4);
	string filename = "tests/test1/liste.txt";
	g->initPatterns("tests/test2/pattern2.txt");
	g->initFromList(filename);

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
//	cout << "Fin fillFreq\n";
//	cout << "Debut Impression\n";
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
	delete g;
}

void doTest1()
{
	cerr << "\n**test dimension tableau**\n";
	FreqKmer *f = new FreqKmer(6);
	string filename = "tests/test3/test.fasta";

	f->initPatterns("tests/test3/3patterns.txt");
	f->initFromFasta(filename);
	if(f->getNCol()!=84 || f->getNLigne()!=107)
	{
		cerr << "Test1 fail...\n";
	}
	else
	{
		cerr << "Test1 OK !\n";
	}

	delete f;
}

void callTest(int i)
{
	switch (i) {
		case 1:
			doTest1();
			break;
		case 2:
			doTest2();
			break;
		case 3:
			doTest3();
			break;
		default:
			break;
	}
}



void executeTests(int k)
{
	cerr << "===Début des tests===\n";
	for(int i=1; i <= k ; i++)
	{
		callTest(i);
	}
	cerr << "\n===Fin des tests===\n";
}
int main(int argc, char **argv) {

	if (doTest)
		executeTests(3);
	//FreqKmer *f = new FreqKmer(50);
//	f->initFromFasta("/home/jeremy/Bureau/boug.fasta");
//	f->initPatterns("/home/jeremy/mitomer/trunk/generate_learn/pattern_test");
	//f->fillFreq();
	//f->imprimeCSV();

	return 0;
}





