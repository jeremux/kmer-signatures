/*
 * count.c
 * 
 * Copyright 2014 Jeremy FONTAINE <jeremux.f@gmail.com>
 * 
 */
 
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "count_lib.h"



int premier = 0;
int taille_premier = 0;
int dernier = 0;
unsigned char dernier_char = 0;
int tmax = 0;


/**
 * Mapping pour transformer un nucletotide 
 * en une valeur
 * pour les caractère a et A => 0
 * 					  c et C => 1
 *                    g et G => 2
 * 					  t et T => 3
 * pour les autres : 4
 * */
const unsigned char map[256] =
{4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4,
4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2,
4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};

const char map_inverse[4] = { 'A', 'C', 'G', 'T' };

/**
 * 
 * @param x 
 * @return 4^x
 * */
int puissance4(int x) 
{
	int t = (int)1 << (x * 2);
	
	if ( t > tmax)
	{
		tmax = t;
	}
	return (int)1 << (x * 2);
}

/**
 * Vérifie qu'un pointeur est non null
 * Si erreur est different de null et ptr est null 
 * @return le message erreur.
 * */ 
void my_error(void *ptr, const char *erreur) 
{
	if (ptr == NULL) 
	{
		if(erreur != NULL)  
		{
			fprintf(stderr, "Erreur: %s - %s\n", erreur, strerror(errno));
		}
		else 
		{
			fprintf(stderr, "Erreur: %s\n", strerror(errno));
		}
		exit(EXIT_FAILURE);
	}
}

/**
 * Retourne val kmer
 * 
 * */
 char *indice_to_kmer(int indice, long kmer)  
 {

	size_t i = 0;
	size_t j = 0;
	char *tab_numerique = calloc(kmer,  sizeof(char));
	char *resultat = calloc(kmer , sizeof(char));
	my_error(tab_numerique, "Erreur calloc indice_to_kmer");
	my_error(resultat, "Erreur calloc indice_to_kmer");
	
	/* Meme principe que le passage du decimal vers binaire */
	while (indice != 0) 
	{
		tab_numerique[i] = indice % 4;
		indice /= 4;
		i++;
	}
	
	for(j = 0; j < (kmer - i); j++)
		resultat[j] = 'A';

	int offset = j;
	
	size_t start = i ;

	i--;
	j = 0;

	// reverse the array, as j increases, decrease i
	for(j = 0; j < start; j++, i--)
		resultat[j + offset] = map_inverse[(int)tab_numerique[i]];
	
  // set our last character to the null termination byte
	resultat[kmer ] = '\0';

	free(tab_numerique);
	return resultat;
}


/**
 *
 *
 *
 **/
 void init(int *tab_comptage,int taille_kmer)
{

	const unsigned long taille = puissance4(taille_kmer);
	int i;

	for(i=0;i<taille+1;i++)
	{
		tab_comptage[i]=0;
	}
}

/**
 *
 *
 *
 *
 **/
 void affiche(int *tab_comptage,int taille_kmer)
{

	const unsigned long taille = puissance4(taille_kmer);
	int i;

	for(i=0;i<taille+1;i++)
	{
		if (tab_comptage[i]!=0)
		{
			printf("%s = %d\n",indice_to_kmer(i,taille_kmer),tab_comptage[i]);
		}
	}
}


void count(const char *seq,const size_t seq_taille,int *tab_bool,int *taille_kmer,int **resultat,int taille_sous_sequence)
{


	int i,k;
	int flag_tour = 0;
	int kmer_taille = 0;
	int pattern_taille = 0;

	int val;

	char *seq_numerique = (char *) malloc(seq_taille*(sizeof(char)));
	my_error(seq_numerique,"Erreur malloc copie seq"); 
	
	/* piste pour stat sur kmer */
	// char *bad_kmer = (char *) malloc(kmer_taille*(sizeof(char)));
	// my_error(bad_kmer,"Erreur malloc string pour kmer inconnu");


	
	/* Pré traitement */
	for(k = 0; k < seq_taille; k++)
	{
		seq_numerique[k] = map[(int)seq[k]];
	}


	val = tab_bool[pattern_taille];

	

	while(val!=-1)
	{
		if (val==1)
		{
			kmer_taille++;
		}
		
		pattern_taille++;

		val = tab_bool[pattern_taille];
	} 
	
	*taille_kmer = kmer_taille;

	// printf("pattern_taille = %d\n",pattern_taille);
		
	
	for(i=0;i<=seq_taille - taille_sous_sequence;i++)
	{


		if (flag_tour)
		{

			// printf("kmer_taille = %d\n",kmer_taille);
			memcpy(resultat[i],resultat[i-1],puissance4(kmer_taille)*sizeof(int));
			// printf("acces à premier = %d\n",premier);
			resultat[i][premier] = resultat[i][premier] - 1;
			add_one(seq_numerique,i,taille_sous_sequence,pattern_taille,resultat[i],tab_bool,kmer_taille);

		}
		else
		{	
			count_seq(seq_numerique,taille_sous_sequence,tab_bool,resultat[i],pattern_taille);
			flag_tour = 1;
		}
	}

	// printf("tmax = %d\n",tmax);
	free(seq_numerique);
}



void add_one(const char *seq,int i,const size_t seq_taille,int pattern_taille,int *resultat,int *tab_bool,int kmer_taille)
{
	int alpha =  i + seq_taille - (pattern_taille + 1);
	int beta = i + seq_taille - 1 ;
	int k;
	int l=0;

	
	// printf("seq[%d] = %d\n",beta,seq[beta]);
	// printf("dernier before = %d\n",dernier);
	if (kmer_taille==pattern_taille)
	{
		dernier = (dernier - (seq[alpha]*puissance4(kmer_taille-1)))*4 + seq[beta];
	}
	else
	{
		dernier = 0;
		for (k = i+seq_taille-pattern_taille,l=0; k < i+seq_taille; k++)
		{

			if(tab_bool[k-(i+seq_taille-pattern_taille)])
			{
				dernier += seq[k] * puissance4(kmer_taille-l-1);
				l++;
				// printf("seq[k] = %d\n",seq[k]);
			}
		}
		// printf("dernier = %d\n",dernier);
	}
	
	// printf("dernier = %d\n",dernier);
	resultat[dernier] = resultat[dernier] + 1;

	premier = 0;
	for (k = 0,l=0; k < pattern_taille; k++)
	{

		if(tab_bool[k])
		{
			premier += seq[i+k] * puissance4(kmer_taille-l-1);
			l++;
		}
	}
	

}



/**
 * Permet le comptage des kmer dans une sequence
 * @param seq: La sequence où compter
 * @param seq_taille: taille de la séquence
 * @param kmer_taille: la taille du kmer
 * @param tab_comptage: la table de comptage
 * */
void count_seq(const char *seq, const size_t seq_taille, int *tab_bool, int *resultat,int pattern_taille) 
{
	long long position;
	long long i;
	
	int lk, lp;
	int cpt = 0;
	int flag_premier = 1;
	int indice_kmer = 0;
	//~ int flag_inconnu = 0;
	
	/* On travaille sur une copie, car en param const char */

			
	/* Parcours de la sequence */
	for(position = 0; position < (signed)(seq_taille - pattern_taille + 1); position++) 
	{
		indice_kmer = 0;

		
		dernier_char = seq[position + pattern_taille - 1];

		/* Pour chaque kmer */
		for(i = position + pattern_taille - 1, lk = 0 ,lp = 0; i >= position; i--, lp++)
		{
			/* si nucléotide inconnu */
			if(seq[i] == 4) 
			{
				/* on va au prochain kmer */
				position = i;
				
				/* goto bad_nuc : on ne compte pas */
				goto bad_nuc;
			}


			/* calcul de l'indice */	
			if(tab_bool[pattern_taille-(lp+1)])
			{	
				// printf("tab_bool[%d] = %d\n", pattern_taille-(lp+1),tab_bool[pattern_taille-(lp+1)]);
				indice_kmer += seq[i] * puissance4(lk);
				lk++;
			}
			
		}
		
		
		/* comptage du kmer */
		cpt++;
		// printf("kmer %d et indice = %d\n",cpt,indice_kmer);
		if(flag_premier)
		{
			premier= indice_kmer;
			taille_premier = (premier == 0 ? 1 : (int)(log10(premier)+1));
			flag_premier = 0;
		}

		// printf("resultat avant = %d\n",resultat[indice_kmer]);
		resultat[indice_kmer]++;
		// printf("resultat apres = %d\n",resultat[indice_kmer]);

		bad_nuc: ;
		

	}
		
	dernier = indice_kmer;

	// printf("=======================\n");
	// printf("Premier = %d\nDernier = %d\nDernier_char = %d\n",premier,dernier,dernier_char);
	// printf("=======================\n");
		// bad_nuc: ;
}

