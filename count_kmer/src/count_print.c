/*
 * count_print.c
 * 
 * Copyright 2014 Jeremy FONTAINE <jeremux.f@gmail.com>
 * 
 */
 
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "count_print.h"
#include "count_lib.h"



void imprime_kmer(FILE *out,int taille_kmer)
{

	int i;

	unsigned long long taille =  puissance4(taille_kmer) ;

	for (i = 0; i < taille; i++)
	{
		fprintf(out,"%s;", indice_to_kmer(i,taille_kmer));
	}

	
}

void imprime_csv(const char *path,int nb_kmer,int *taille_kmer,int nombre_sous_sequence,int ***resultat)
{
	int i,j,k;
	FILE *out = NULL;
	unsigned long long nb_possibilite;


	out = fopen(path, "a");
	if (out != NULL)
	{
		for(i=0;i<nb_kmer;i++)
		{
			imprime_kmer(out,taille_kmer[i]);
		}
		
		fprintf(out, "\n");

		for(i=0;i<nombre_sous_sequence;i++)
		{
			// printf("J'itère\n");
			for( k= 0;k<nb_kmer;k++)
			{
				nb_possibilite = puissance4(taille_kmer[k]);
				// printf("nb_possibilite = %llu\n",nb_possibilite);

				for (j = 0; j < nb_possibilite; j++)
				{
					fprintf(out,"%d;",resultat[k][i][j]);
				}
				
			}
			fprintf(out, "\n");
		}

		fclose(out);
	}
	else
	{
		fprintf(stderr, "Erreur creation du fichier %s\n",path);
		exit(EXIT_FAILURE);
	}
		
}

void imprime_entete_weka(FILE *fi,int taille_kmer)
{
	int i;

	unsigned long long taille =  puissance4(taille_kmer) ;

	for (i = 0; i < taille; i++)
	{
		fprintf(fi,"@attribute %s numeric\n", indice_to_kmer(i,taille_kmer));
	}
}

void imprime_weka(const char *path,int nb_kmer,int *taille_kmer,int nombre_sous_sequence,int ***resultat,const char *taxid,int flag_en_tete,char *les_accessions)
{
	int i,j,k;
	FILE *out = NULL;
	unsigned long long nb_possibilite;

	// printf("JECRIS\n");
	out = fopen(path, "a");
	if (out != NULL)
	{	
		// if (flag_en_tete==0)
		// {
		// 	fprintf(out, "@relation kmers_count\n");
		// 	for(i=0;i<nb_kmer;i++)
		// 	{
		// 		imprime_entete_weka(out,taille_kmer[i]);
		// 	}

		// 	fprintf(out,"@attribute id { }\n");
		// 	fprintf(out, "@data\n");
		// }


		for(i=0;i<nombre_sous_sequence;i++)
		{
			// printf("J'itère !");
			for( k= 0;k<nb_kmer;k++)
			{
				nb_possibilite = puissance4(taille_kmer[k]);
				// printf("nb_possibilite = %llu\n",nb_possibilite);

				for (j = 0; j < nb_possibilite; j++)
				{
					fprintf(out,"%d,",resultat[k][i][j]);
					// fprintf(stdout,"%d,",resultat[k][i][j]);

				}
				
			}
			fprintf(out,"%s %% %s\n",taxid,les_accessions);
			// fprintf(stdout,"%s %% %s\n",taxid,les_accessions);
		}

		fclose(out);
	}
	else
	{
		fprintf(stderr, "Erreur creation du fichier %s\n",path);
		exit(EXIT_FAILURE);
	}	
}

