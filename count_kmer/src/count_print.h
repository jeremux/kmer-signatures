#ifndef COUNT_PRINT_H_
#define COUNT_PRINT_H_


void imprime(int **tab_comptage,int *taille_kmer);

void imprime_weka(const char *path,int nb_kmer,int *taille_kmer,int nombre_sous_sequence,int ***resultat,const char *taxid,int flag,char *les_accessions);

void imprime_csv(const char *path,int nb_kmer,int *taille_kmer,int nombre_sous_sequence,int ***resultat);

void imprime_kmer(FILE *out,int taille_kmer);

void imprime_entete_weka(FILE *fi,int taille_kmer);

char *indice_to_kmer(int indice, long kmer);

#endif