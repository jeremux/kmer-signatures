/*
 * classPattern.cpp
 *
 *  Created on: 10 juil. 2014
 *      Author: jeremy
 */

#include "classPattern.h"
#include <algorithm>
#include <cmath>

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

string Pattern::getKmer(string seq,int coord)
{
	string res = "";
	for (int var = 0; var < getTaillePattern() ; var++)
	{
		if(extraire(var))
		{
			res += seq[coord];
		}
		coord++;
	}
	return res;
}


Pattern::~Pattern() {

}

