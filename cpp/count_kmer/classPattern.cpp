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


Pattern::~Pattern() {

}

