//============================================================================
// Name        : kmer_frequencies.cpp
// Author      : Jeremy FONTAINE
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <istream>
#include <sstream>

using namespace std;

int main() {
	string s = "hello world";
	istringstream iss(s);
	string token1, token2;
	iss >> token1;
	iss >> token2;
	cout << token1 << " " << token2 << endl; // prints
	return 0;
}
