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

void doTest2()
{
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
		cerr << "Test2 OK !\n";
	}
	else
	{
		cerr << "Test2 fail...\n";
	}

}

void doTest1()
{
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
		cerr << "Test1 OK !\n";
	}
	else
	{
		cerr << "Test1 fail...\n";
	}

}


void executeTests()
{
	cerr << "===Début des tests===\n";
	cerr << "**test pattern simple**\n";
	doTest1();
	cerr << "**test multi-pattern**\n";
	doTest2();
	cerr << "===Début des tests===\n";
}
int main(int argc, char **argv) {

	if (doTest)
		executeTests();

//	for(int j=0;j<col;j++)
//		{
//			cout << j << ";";
//		}
//		cout << "\n";
//		for(int i=0;i<ligne;i++)
//		{
//			for(int j=0;j<col;j++)
//			{
//				cout << freq[i][j] << ";";
//			}
//			cout << "\n";
//		}
//	f->imprimeCSV();
//	cout << "Fin Impression\n";
	return 0;
}





