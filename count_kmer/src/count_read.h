#ifndef COUNT_READ_H_
#define COUNT_READ_H_



int **make_tab_bool(int nb_kmer,int taille_max_kmer);

int read_pattern(const char *path,int **tab_bool);

// char *read_seq(const char *path, int *nb,char **seq);

char *read_seq(const char *path, int *nb,char **seq,char **acc);


void free_tab_bool(int **tab,int nb_kmer);

#endif