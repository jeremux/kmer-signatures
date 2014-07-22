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

void doTest10()
{
        cerr << "\n**test decalage bis**\n";

        FreqKmer *f = new FreqKmer(3,2);
        string filename = "tests/test1/liste2.txt";
        f->initPatterns("tests/test1/pattern2.txt");
        f->initDataFromListFastaPath(filename);



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

        if (f->getNLigne()!=2 || f->getNCol()!=16)
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
                		cerr << "freq[" << i << "][" << j << "] = " << f->getFreq()[i][j] <<  " || fb[" << i << "][" << j << "] = " << freq[i][j]<< "\n";
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
        delete f;
}

void doTest9()
{
	cerr << "\n**test decalage 2**\n";

	        FreqKmer *f = new FreqKmer(5,2);
	        string filename = "tests/test5/liste.txt";
	        f->initPatterns("tests/test5/pattern.txt");
	        f->initDataFromListFastaPath(filename);



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


	        f->fillFreq();
	        f->imprimeCSV("out.csv");

	        if (f->getNLigne()!=5 || f->getNCol()!=20)
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
	                        if(f->getFreq()[i][j]!=freq[i][j])
	                        {
	                                res=false;
	                                break;
	                        }
	                }
	        }
	        printResult(res,9);
	        delete f;
}

void doTest8()
{
        cerr << "\n**test decalage**\n";

        FreqKmer *f = new FreqKmer(3,2);
        string filename = "tests/test1/liste.txt";
        f->initPatterns("tests/test1/pattern.txt");
        f->initDataFromListFastaPath(filename);



        int ligne = 5;
        int col = 20;
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

        if (f->getNLigne()!=5 || f->getNCol()!=20)
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
                		cerr << "freq[" << i << "][" << j << "] = " << f->getFreq()[i][j] <<  " || fb[" << i << "][" << j << "] = " << freq[i][j]<< "\n";
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
        delete f;
}

void doTest7()
{
	cerr << "\n**Test de la map deux**\n";
	FreqKmer *f = new FreqKmer(4);
	bool res = true;
	string filename = "tests/test4/liste.txt";
	f->initPatterns("tests/test4/pattern.txt");
	f->initDataFromListFastaPath(filename);

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
}

void doTest6()
{
	cerr << "\n**Test de la map une**\n";
	FreqKmer *f = new FreqKmer(4);
	bool res = true;
	string filename = "tests/test1/liste.txt";
	f->initPatterns("tests/test1/pattern.txt");
	f->initDataFromListFastaPath(filename);

	res = !(f->obtainNbLineData(0)!=3 || f->obtainNbLineData(1)!=4 || f->obtainStartLineData(0)!=0 || f->obtainEndLineData(0)!=2 || f->obtainStartLineData(1)!=3 || f->obtainEndLineData(1)!=6);

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
	f->initDataFromListFastaPath(filename);

	if (f->obtainStartColKmer(0)!=0 || f->obtainStartColKmer(1)!=16 || f->obtainEndColKmer(0)!=15 || f->obtainEndColKmer(1)!=19)
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
	f->initDataFromListFastaPath(filename);

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
        f->initDataFromListFastaPath(filename);

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
        g->initDataFromListFastaPath(filename);

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
				case 8:
                		doTest8();
                		break;
				case 9:
                		doTest9();
                		break;
				case 10:
						doTest10();
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


