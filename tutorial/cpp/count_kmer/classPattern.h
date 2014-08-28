/*!
 * \file classPattern.h
 * \brief to manage kmer
 * \author jeremy FONTAINE
 * \version 1.0
 */

#ifndef CLASSPATTERN_H_
#define CLASSPATTERN_H_


#include <string>
#include <vector>
using namespace std;

/*! \class Pattern
* \brief K-mer object
*
*/
class Pattern {

private:
	string pattern; /*!< the pattern: ex ##_## */
	int patternSize; /*!< pattern length */
	int kmerSize; /*!< kmer length */

public:

	/*!
	 *  \brief Constructeur
	 *  \param p : le pattern
	 */
	Pattern(string p);

	/*!
	 *  \brief Destructeur
	 *
	 *  Destructeur de la classe Pattern
	 */
	virtual ~Pattern();

	/*!
	 *  \brief Le nombre de kmer
	 *
	 *  Methode qui permet d'avoir le nombre de
	 *  kmer possible selon la taille du kmer
	 *
	 *  \return le nombre de kmer du pattern courant
	 */
	int getAllCombi();


	int getSizePattern(){return patternSize;}

	int getSizeKmer(){return kmerSize;}

	/**
	 * Permet d'avoir la liste des kmer
	 * @return liste des kmer: a...a,a...c,...,t..tg,t..tt
	 */
	vector<string> getCombi();

	/**
	 * @return vrai si c'est un pattern continue
	 */
	bool isContinue();


	int getKmer(int* seq,int coord);

	/**
	 * Permet de savoir si on
	 * doit considerer une position
	 * dans le pattern pour le
	 * comptage
	 * @param i: position dans le kmer
	 */
	bool extraire(int i);
};

#endif /* CLASSPATTERN_H_ */
