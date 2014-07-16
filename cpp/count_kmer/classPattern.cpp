/*
 * classPattern.cpp
 *
 *  Created on: 10 juil. 2014
 *      Author: jeremy
 */

#include "classPattern.h"
#include <algorithm>
#include <cmath>
#include<iostream>
#include <vector>

using namespace std;

Pattern::Pattern() {
	// TODO Auto-generated constructor stub
	this->pattern = "";
}

Pattern::Pattern(string p)
{
	this->pattern = p;
}

int Pattern::getTaillePattern()
{
	int res;
	res = this->pattern.size();
	return res;
}


int Pattern::getTailleKmer()
{
	int res;
	res = count(this->pattern.begin(),pattern.end(),'#');
	return res;
}

int Pattern::getAllCombi()
{
	int res;
	res = pow(4,getTailleKmer());
	return res;
}

bool Pattern::isContinue()
{
	return (getTailleKmer()==getTaillePattern());
}

bool Pattern::extraire(int i)
{
	return (this->pattern[i]=='#');
}

int Pattern::getKmer(int* seq,int coord)
{
	int res = 0;
	int cpt = 0;
	for (int var = 0; var < getTaillePattern() ; var++)
	{
		if(extraire(var))
		{
			res += seq[coord]*pow(4,getTailleKmer()-1-cpt);
			cpt++;
		}
		coord++;
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
	for(int i=0;i<getTailleKmer()-1;i++)
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

}

