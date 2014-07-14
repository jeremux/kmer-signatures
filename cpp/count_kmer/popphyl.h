#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <typeinfo>
#include <iomanip>

#include <pwd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <vector>
#include <libgen.h>




// Deprecated
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <time.h>


using namespace std;

#ifndef CONSTANTS
#define CONSTANTS

const double PI = 3.14159265358979;
const int TAILLE_FENETRE=100;

const int nAA=20;
const int     MAXSITE=1000000;
const char    AminoAcids[] = {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','-', 'X'};
const string  AAcolors[] = {"blue", "yellow", "magenta", "magenta", "cyan", "black", "cyan", "blue", "red", "blue", "blue", "green", "yellow", "green", "red", "green", "green", "blue", "cyan", "cyan", "black"};

//                            A , C , D , E , F , G , H , I , K , L , M , N , P , Q , R , S , T , V , W , Y
const int taylorcat[][20] ={{ 1 , 1 , 0 , 0 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 0 , 0 , 0 , 0 , 0 , 1 , 1 , 1 , 1 }, // hydrophobic
			    			{ 0 , 0 , 0 , 0 , 0 , 0 , 1 , 0 , 1 , 0 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 , 0 }, // positive
			    			{ 0 , 0 , 1 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }, // negative
			    			{ 0 , 0 , 1 , 1 , 0 , 0 , 1 , 0 , 1 , 0 , 0 , 1 , 0 , 1 , 1 , 1 , 1 , 0 , 1 , 1 }, // polar
			    			{ 0 , 0 , 1 , 1 , 0 , 0 , 1 , 0 , 1 , 0 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 , 0 }, // charged
			    			{ 1 , 1 , 1 , 0 , 0 , 1 , 0 , 0 , 0 , 0 , 0 , 1 , 1 , 0 , 0 , 1 , 1 , 1 , 0 , 0 }, // small
			    			{ 1 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 }, // tiny
			    			{ 0 , 0 , 0 , 0 , 1 , 0 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 1 }, // aromatic
			    			{ 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 0 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 0 , 0 }};// aliphatic
//                   sum      3   2   4   3   2   3   5   2   4   2   1   2   1   1   3   3   3   3   3   3
const string taylorname[]   ={"hydrophobic", "positive", "negative", "polar", "charged", "small", "tiny", "aromatic", "aliphatic"};
const int taylorn=9;
const int taylorfreedom[]   ={13,3,2,11,5,9,3,4,3};




const int ncolors=7;
const string  colors[] = {"red", "blue", "magenta", "green", "cyan", "yellow", "black"};

const int     nNT=4;
const char    DNAletters[] = {'A','C','G','T','-', 'N'};
const char    dnaletters[] = {'a','c','g','t','-', 'n'};
const string  NTcolors[] = {"red", "green", "yellow", "blue", "black", "black"};
const double  NTpstrickScale[] = {1, 1, 1, 1.065};  // letter scale correction for vertical logo
const double  NTpstrickShift[] = {0, 0.01, 0.01, 0};      // letter position correction for vertical logo

const int NautorizedGeneName=13;
const char* const autorizedGeneName[NautorizedGeneName]={"CYTB", "COX1", "COX2", "COX3", "ND1", "ND2", "ND3", "ND4", "ND5", "ND6", "ND4L", "ATP8", "ATP6"};


enum datatype{Nexus=0, Phylip=1, Fasta=2};
enum Switch{Yes=1,No=0};
enum Alphabet{dna=0,protein=1};




static string printSwitch(Switch s)
{
  if(s==Yes)
    return "Yes";
  else
    return "No";
}

static Switch getSwitch(int s)
{
  if(s==1)
    return Yes;
  else
    return No;
}

const bool dataVerbose=false;

#endif
