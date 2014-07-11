/*
 * classPattern.h
 *
 *  Created on: 10 juil. 2014
 *      Author: jeremy
 */

#ifndef CLASSPATTERN_H_
#define CLASSPATTERN_H_
#include <string>
using namespace std;

class Pattern {

private:
	string pattern;

public:
	Pattern();
	Pattern(string p);
	virtual ~Pattern();

	int getAllCombi();
	int getTaillePattern();
	int getTailleKmer();
	bool isContinue();
	string getKmer(string seq,int coord); /* coord de 0 à taille(seq)-1 */
	bool extraire(int i);
};

#endif /* CLASSPATTERN_H_ */
