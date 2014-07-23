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

void doTest11()
{
        cerr << "\n**Données réelles validées par perl/C**\n";

        FreqKmer *f = new FreqKmer(-1,1);
        string filename = "tests/haemo/haemo.fasta";
        f->initPatterns("tests/haemo/pattern5.txt");
        f->initFromFasta(filename);

        int val[] = {11,10,7,18,12,5,5,7,10,3,7,6,29,4,14,17,8,7,9,14,7,2,3,6,4,2,3,5,4,6,8,13,14,4,9,7,2,3,1,5,5,1,
        		0,9,7,1,3,9,25,14,8,31,4,8,1,4,10,8,6,19,24,8,13,21,8,5,7,10,9,3,1,4,7,5,7,4,16,5,7,19,9,4,4,11,
        		4,2,1,1,1,3,1,9,7,4,3,10,5,1,4,3,1,2,0,4,0,0,0,3,6,1,1,4,10,7,4,12,8,1,2,6,1,3,6,6,3,8,5,13,9,7,7,
        		11,2,2,4,2,6,1,2,3,13,4,10,6,5,3,6,6,7,1,3,6,1,0,1,2,7,5,4,4,8,4,5,10,2,1,2,1,2,1,0,1,5,1,2,8,5,2,7,
        		19,2,1,2,1,3,2,4,4,11,7,5,9,28,8,3,33,16,11,3,15,10,10,6,11,34,11,16,46,8,3,6,13,7,3,1,10,2,4,1,1,1,
        		5,6,7,6,4,7,15,13,6,1,14,20,2,8,7,21,8,5,20,27,23,17,38,15,8,2,22,8,8,8,7,30,19,18,25,9,3,8,12,7,1,1,
        		10,13,5,3,3,11,4,9,13,4,4,0,8,3,0,3,4,2,1,0,2,5,2,1,4,8,3,0,9,5,6,0,5,8,4,2,3,9,2,3,3,7,9,8,20,6,4,3,
        		4,4,10,8,9,16,14,3,19,8,4,8,8,3,3,0,4,5,3,5,5,6,3,9,13,5,0,5,3,2,0,0,1,1,0,0,0,0,0,2,3,0,2,0,0,0,2,0,2,
        		1,0,2,1,5,3,0,4,2,4,1,12,3,2,1,0,7,0,4,4,3,2,4,12,0,4,0,4,3,0,0,3,0,3,0,2,2,1,1,2,1,2,2,1,1,1,2,2,0,0,0,
        		1,1,1,3,4,2,0,1,0,0,0,0,1,0,2,1,0,3,1,1,3,5,3,2,5,2,1,1,2,1,0,0,1,6,1,1,7,12,5,1,6,8,7,1,5,1,5,2,2,15,1,
        		9,15,7,3,2,12,4,1,0,1,2,1,1,3,7,5,0,4,6,2,2,6,1,0,1,2,3,2,3,8,7,6,3,6,5,5,2,11,7,7,2,6,2,0,2,6,16,6,8,
        		15,8,6,1,7,7,1,2,9,3,1,2,8,14,3,2,9,5,0,1,3,1,1,4,2,1,4,0,2,4,2,3,4,2,3,1,5,3,3,0,3,4,0,0,1,5,2,3,8,11,
        		8,6,16,3,3,1,3,4,7,5,11,14,5,3,11,5,4,5,6,2,2,0,2,5,3,1,3,7,3,3,6,4,0,3,7,2,0,0,1,0,1,3,2,6,1,2,3,1,1,0,
        		1,0,0,0,0,0,0,0,2,1,0,1,3,4,4,1,7,6,2,1,3,2,0,3,5,7,6,0,9,6,4,3,7,2,1,1,5,3,4,0,6,14,1,9,16,2,0,1,2,2,1,
        		0,1,0,0,1,2,3,1,1,1,3,3,1,4,1,2,1,1,3,1,1,1,1,2,1,4,7,2,2,9,3,1,2,3,1,4,4,3,6,7,4,7,9,6,3,11,9,3,0,6,8,0,
        		1,8,13,7,9,22,4,1,3,6,4,0,0,1,1,1,2,3,1,3,1,4,2,3,2,2,3,1,0,3,5,0,1,5,6,0,
        		0,6,9,5,3,18,13,5,3,6,6,3,4,3,16,7,3,7,18,11,10,27,12,11,7,7,8,2,3,3,26,4,18,26,13,6,
        		13,23,17,5,4,12,7,0,0,3,21,7,4,8,10,0,2,12,10,5,3,7,10,1,2,3,12,1,3,12,29,14,13,40,17,7,2,8,14,8,18,
        		15,51,20,12,40,11,6,5,14,1,2,4,3,2,5,4,7,16,7,11,14,10,6,7,10,5,1,0,2,0,0,0,1,6,1,9,6,1,2,1,2,5,2,1,3,
        		2,1,1,2,4,2,0,4,9,6,4,11,7,1,3,7,6,0,4,6,10,6,2,11,7,4,4,6,2,5,2,3,2,1,3,7,12,4,6,9,12,1,3,10,4,0,1,3,2,
        		0,0,0,5,5,3,12,7,2,6,26,2,1,0,3,6,1,4,6,11,5,7,10,11,10,6,18,7,2,2,4,3,1,3,4,12,12,6,11,16,17,9,24,21,
        		17,6,14,5,10,7,7,34,14,21,39,16,3,7,17,18,4,0,10,1,5,2,3,20,5,9,14,7,3,2,8,9,1,0,6,13,2,5,13,12,0,3,9,24,
        		25,7,40,9,11,4,14,4,4,19,8,34,6,7,14
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
			   delete freq[i];
	   }
        delete freq;

    //    delete f;
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
        printResult(res,10);

        for(int i=0 ; i<ligne ; i++)
	   {
			   delete freq[i];
	   }
        delete freq;

        //delete f;
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

	        if (f->getNLine()!=5 || f->getNCol()!=20)
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

        if (f->getNLine()!=5 || f->getNCol()!=20)
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
        printResult(res,8);
//        delete f;
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
        if(f->getNCol()!=84 || f->getNLine()!=107)
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
				case 11:
						doTest11();
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


