/**
 * \file classPattern.cpp
 * \brief classe permettant de gérer un pattern de kmer.
 * \author Jeremy FONTAINE
 * \version 0.1
 * \date 13/07/2014
 *
 * Cette classe est utilisé pour réprenster
 * un pattern de kmer sous la forme de ##_##
 *
 */

#include "classPattern.h"
#include <algorithm>
#include <cmath>
#include<iostream>
#include <vector>

using namespace std;


Pattern::Pattern(string p)
{
	this->pattern = p;
	this->kmerSize = count(this->pattern.begin(),pattern.end(),'#');
	this->patternSize = this->pattern.size();
}


int Pattern::getAllCombi()
{
	int res;
	res = pow(4,getSizeKmer());
	return res;
}

bool Pattern::isContinue()
{
	return (getSizeKmer()==getSizePattern());
}

bool Pattern::extraire(int i)
{
	return (this->pattern[i]=='#');
}

int Pattern::getKmer(int* seq,int coord)
{
	int res = 0;
	int cpt = 0;

	for (int var = 0; var < getSizePattern() ; var++)
	{
		if(extraire(var))
		{

			res += seq[coord+var]*pow(4,getSizeKmer()-1-cpt);
			cpt++;
		}

	}

	return res;
}

vector<string> Pattern::getCombi()
{
	vector<string> bases;
	bases.push_back("A");
	bases.push_back("C");
	bases.push_back("G");
	bases.push_back("T");
	vector<string> words;
	words = bases;

	vector<string> newords;
	for(int i=0;i<getSizeKmer()-1;i++)
	{
		newords.clear();
		for(size_t j=0;j<words.size();j++)
		{
			for(size_t k=0;k<bases.size();k++)
			{
				newords.push_back(words[j]+bases[k]);
			}
		}
		words.clear();
		words = newords;
	}

//	for(size_t l=0;l<words.size();l++)
//	{
//		cerr << "@attribute " << words[l] << " numeric \n";
//	}

	return words;
}
Pattern::~Pattern() {
	pattern="";
	kmerSize=0;
	patternSize=0;

}

