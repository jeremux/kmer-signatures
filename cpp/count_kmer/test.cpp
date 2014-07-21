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

void doTest7()
{
	cerr << "\n**Test de la map deux**\n";
	FreqKmer *f = new FreqKmer(4);
	bool res = true;
	string filename = "tests/test4/liste.txt";
	f->initPatterns("tests/test4/pattern.txt");
	f->initFromList(filename);

	res = (
		   f->getStartLineDataSeq(0,0)==0
		&& f->getEndLineDataSeq(0,0)==2
		&& f->getNbLineDataSeq(0,0)==3

		&& f->getStartLineDataSeq(0,1)==3
		&& f->getEndLineDataSeq(0,1)==3
		&& f->getNbLineDataSeq(0,1)==1

		&& f->getStartLineDataSeq(1,0)==4
		&& f->getEndLineDataSeq(1,0)==7
		&& f->getNbLineDataSeq(1,0)==4

		&& f->getStartLineDataSeq(1,1)==8
		&& f->getEndLineDataSeq(1,1)==9
		&& f->getNbLineDataSeq(1,1)==2

		&& f->getStartLineDataSeq(1,2)==10
		&& f->getEndLineDataSeq(1,2)==14
		&& f->getNbLineDataSeq(1,2)==5
	);



	printResult(res,7);

	delete f;
}

void doTest6()
{
	cerr << "\n**Test de la map une**\n";
	FreqKmer *f = new FreqKmer(4);
	bool res = true;
	string filename = "tests/test1/liste.txt";
	f->initPatterns("tests/test1/pattern.txt");
	f->initFromList(filename);

	res = !(f->getNbLineData(0)!=3 || f->getNbLineData(1)!=4 || f->getStartLineData(0)!=0 || f->getEndLineData(0)!=2 || f->getStartLineData(1)!=3 || f->getEndLineData(1)!=6);

	printResult(res,6);

	delete f;
}
void doTest5()
{
	cerr << "\n**Test espace des kmers**\n";
	FreqKmer *f = new FreqKmer(4);
	bool res = true;
	string filename = "tests/test1/liste.txt";
	f->initPatterns("tests/test1/pattern.txt");
	f->initFromList(filename);

	if (f->getStartColKmer(0)!=0 || f->getStartColKmer(1)!=16 || f->getEndColKmer(0)!=15 || f->getEndColKmer(1)!=19)
	{
		res=false;
	}

	printResult(res,5);

	delete f;
}
void doTest4()
{
	cerr << "\n**Test taille des sequences**\n";
	FreqKmer *f = new FreqKmer(4);
	bool res = true;
	string filename = "tests/test1/liste.txt";
	f->initPatterns("tests/test1/pattern.txt");
	f->initFromList(filename);

	if (f->getData()[0]->getLengthSeq(0)!=6 || f->getData()[1]->getLengthSeq(0)!=7)
	{
		res=false;
	}

	printResult(res,4);

	delete f;
}
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
//      cout << "Fin fillFreq\n";
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
                case 4:
                		doTest4();
                		break;
                case 5:
                		doTest5();
                		break;
                case 6:
                		doTest6();
                		break;
                case 7:
                		doTest7();
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


