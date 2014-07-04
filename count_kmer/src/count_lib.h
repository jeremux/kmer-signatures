#ifndef COUNT_LIB_H_
#define COUNT_LIB_H_


unsigned long long puissance4(unsigned long long x);

/**
 * Vérifie qu'un pointeur est non null
 * Si erreur est different de null et ptr est null 
 * @return le message erreur.
 * */ 
void my_error(void *ptr, const char *erreur);


/**
 * Retourne val kmer
 * 
 * */
 char *indice_to_kmer(unsigned long long indice, long kmer);  

/**
 *
 *
 *
 **/
 void init(unsigned long long *tab_comptage,int taille_kmer);
/**
 *
 *
 *
 *
 **/

void concat_one(unsigned long long *tab_comptage,int taille_kmer,char *dest, int i);

void concat_two(unsigned long long *tab_comptage,char *dest,int j);

void affiche(unsigned long long *tab_comptage,int taille_kmer);

void add_one(const char *seq,int i,const size_t seq_taille,int pattern_taille,int *resultat,int *tab_bool,int kmer_taille);

void count_seq(const char *seq, const size_t seq_taille,int *tab_bool,int *resultat,int pattern_taille) ;

void count(const char *seq,const size_t seq_taille,int *tab_bool,int *taille_kmer,int **resultat,int taille_sous_sequence);

// void count(const char *seq,const size_t seq_taille,const char *kmer_string,int kmer_taille,int pattern_taille,unsigned long long *tab_comptage,int taille_sous_sequence);



#endif 

