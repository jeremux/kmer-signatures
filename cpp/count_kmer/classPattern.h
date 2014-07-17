/*!
 * \file classPattern.h
 * \brief classe pour gérer les kmers
 * \author jeremy FONTAINE
 * \version 1.0
 */

#ifndef CLASSPATTERN_H_
#define CLASSPATTERN_H_
#include <string>
#include <vector>
using namespace std;

/*! \class CPlayer
* \brief classe representant un pattern de kmer
*
*  La classe gere la taille, les combinaisons...
*/
class Pattern {

private:
	string pattern; /*!< Le pattern */

public:
	Pattern();
	/*!
	 *  \brief Constructeur
	 *
	 *  Constructeur de la classe Pattern
	 *
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
	int getTaillePattern();
	int getTailleKmer();
	vector<string> getCombi();
	bool isContinue();
	int getKmer(int* seq,int coord); /* coord de 0 à taille(seq)-1 */
	bool extraire(int i);
};

#endif /* CLASSPATTERN_H_ */
