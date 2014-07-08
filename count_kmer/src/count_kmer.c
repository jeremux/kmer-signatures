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
#include <getopt.h>
#include "count_lib.h"
#include "count_read.h"
#include "count_print.h"


#define TAILLE_MAX_SEQ 100000 /* TAILLE MAX DUNE SEQUENCE */
#define MAX_SEQ 20000         /* NOMBRE DE SEQUENCES MAX */
#define VERSION 1.0
#define TAILLE_FENETRE 100    /* TAILLE PAR DEFAUT DES READS */
#define TAILLE_STRING_ACCESSION 20


void help() {
	printf("usage: count_kmer -i fasta_file -k pattern_kmer -l x -o file_out \n\n"
				 "count mers in fasta_file with specified kmer in file pattern_kmer\n"
				 "\n"
				 "  --input    -i  fasta file with sequences to count\n"
				 "  --kmer     -k  file with pattern of mers to count\n"
				 "  --length   -l  length of the read\n"
				 "  --output   -o  output filename\n"
				 "  --taxid    -t  taxid of the taxon\n"
				 "  --version  -v  version of the program\n"
				 "  --help     -h  print this help\n"
				 "\n"
				 " jeremy.fontaine@etudiant.univ-lille1.fr\n");
}


int main(int argc,char *argv[])
{
	
	int i,j,val,x,nb_kmer;
	// FILE *tmp;
	int taille_fenetre = 0;
	int nb_sequences;
	char *path_fasta = NULL;
	char *path_kmer = NULL;
	char *out = NULL;
	char *taxid = NULL;
	int ***resultat = NULL;
	// char *taxid = NULL;
	int flag_fenetre = 1;
	int taille_fenetre_origine;
	// int flag_en_tete = 1;
	

	int nb_possibilite;
	int nombre_sous_sequence ;

	

	static struct option long_options[] = {
		{"input", required_argument, 0, 'i'},
		{"kmer",  required_argument, 0, 'k'},
		{"length", required_argument, 0, 'l'},
		{"taxid", required_argument, 0, 't'},
		{"output", required_argument, 0, 'o'},
		{"version", no_argument, 0, 'v'},
		{"help", no_argument, 0, 'h'},
		{0, 0, 0, 0}
	};

	while (1) {


		int option_index = 0;
		int c = 0;

		c = getopt_long (argc, argv, "i:k:l:t:o:vh", long_options, &option_index);

		if (c == -1)
			break;

		switch (c) {
			case 'i':
				path_fasta = optarg;
				break;
			case 'k':
				path_kmer = optarg;
				break;
			case 'l':
				taille_fenetre_origine = atoi(optarg);
				flag_fenetre = 0;
				break;
			case 'o':
				out = optarg;
				break;
			case 't':
				taxid = optarg;
				break;
			case 'h':
				help();
				exit(EXIT_SUCCESS);
				break;
			case 'v':
				printf("count_kmer version %.1f\n",VERSION );
				exit(EXIT_SUCCESS);
				break;
			default:
				break;
		}
	}

	if(argc == 1) {
		help();
		exit(EXIT_FAILURE);
	}

	if(path_fasta == NULL) {
		fprintf(stderr, "Missing input fasta file: count_kmer [-i|--input] seq.fa\n");
		exit(EXIT_FAILURE);
	}

	if(path_kmer == NULL) {
		fprintf(stderr, "Missing input pattern file for mers: count_kmer [-k|--kmer] pattern_kmer\n");
		exit(EXIT_FAILURE);
	}

	if(out == NULL) {
		fprintf(stderr, "Missing output filename\n");
		exit(EXIT_FAILURE);
	}

	if (taxid==NULL)
	{
		fprintf(stderr, "Missing taxid: count_kmer -t the_taxid\n");
		exit(EXIT_FAILURE);
	}
	/*********************************************************************************
	**********************************************************************************
	*******************************PREPROCESS KMER************************************
	**********************************************************************************/

	int **tab_bool = make_tab_bool(15,15); /* macro */

	nb_kmer = read_pattern(path_kmer,tab_bool);

	int *taille_kmer = calloc(nb_kmer,sizeof(int));
	my_error(taille_kmer,"Erreur malloc taille_kmer");

	int *taille_pattern = calloc(nb_kmer,sizeof(int));
	my_error(taille_pattern,"Erreur calloc taille_pattern");
	


	for(i=0;i<nb_kmer;i++)
	{
		j=0;
		int kmer_taille = 0;
		val = tab_bool[i][j];


		while(val!=-1)
		{
			if (val==1)
			{
				kmer_taille++;
			}
			j++;
			val = tab_bool[i][j];
		} 
		
		taille_pattern[i] = j;
		taille_kmer[i] = kmer_taille;
	}
	/*************************************************************************************
	**************************************************************************************
	*******************************PREPROCESS KMER FIN************************************
	**************************************************************************************/


	/*********************************************************************************
	**********************************************************************************
	*******************************READ SEQ*******************************************
	**********************************************************************************/
	char **les_sequences = (char **)(calloc(MAX_SEQ,sizeof(char *)));
	char **les_accessions = (char **)(malloc(MAX_SEQ*sizeof(char *)));


	/* INIT LES SEQ */
	if(les_sequences==NULL)
	{
		fprintf(stderr, "Erreur read_seq (premier malloc)\n");
		exit(EXIT_FAILURE);
	}

	for (i = 0; i < MAX_SEQ; i++)
	{
		// printf("i = %d\n",i);
		printf("Je vais allouer les_sequences[%d]\n",i);
		les_sequences[i] = (char *)(calloc(TAILLE_MAX_SEQ,sizeof(char)));
		if(les_sequences[i]==NULL)
		{
			fprintf(stderr, "Erreur read_seq (second malloc)\n");
			exit(EXIT_FAILURE);
		}
		printf("Fin allocation\n");
		les_sequences[i][0] = '\0';
	}


	/* INIT LES ACC */
	if(les_accessions==NULL)
	{
		fprintf(stderr, "Erreur alloc les_accessions (premier malloc)\n");
		exit(EXIT_FAILURE);
	}


	/* taille max d'une accesion: max */
	for (i = 0; i < MAX_SEQ; i++)
	{
		les_accessions[i] = (char *)(calloc(TAILLE_STRING_ACCESSION ,sizeof(char)));
		if(les_accessions[i]==NULL)
		{
			fprintf(stderr, "Erreur alloc les_accessions (second malloc)\n");
			exit(EXIT_FAILURE);
		}
		les_accessions[i][0] = '\0';
	}

	read_seq(path_fasta,&nb_sequences,les_sequences,les_accessions);



	/*********************************************************************************
	**********************************************************************************
	*******************************READ SEQ FIN***************************************
	**********************************************************************************/

	/*********************************************************************************
	**********************************************************************************
	*******************************ALLOC RESULT **************************************
	**********************************************************************************/

	/*********************************************************************************
	**********************************************************************************
	******************************ALLOC RESULT FIN************************************
	**********************************************************************************/
	/* façon degeulasse de vider le fichier : */

	if (fopen(out,"r")!=NULL)
	{
		// printf("le fichier existe, mode ajout%s\n",out );
	}
	// tmp = fopen(out,"w");
	// fclose(tmp);
	// int indice_ligne = 0;

	/****************************/

	// printf("nb_sequences = %d\n",nb_sequences);
	for(x=0;x<nb_sequences;x++)
	{
		// printf("\n\n***********************************************\n");
		// printf("Traitement sequence %d / %d\n",x+1,nb_sequences);
		// printf("***********************************************\n");
		taille_fenetre = taille_fenetre_origine;

		if(x>0)
		{
			for(i=0; i < nb_kmer; i++)
			{
				for(j=0; j < nombre_sous_sequence; j++)
				{
					free(resultat[i][j]);
				}
				free(resultat[i]);
			}

			free(resultat);
		}

		


		// printf("seq = %s\n",les_sequences[x]);

		int taille_seq = strlen(les_sequences[x]);

		// printf("taille_seq = %d\n",taille_seq);
			

		

		if (flag_fenetre)
		{
			fprintf(stderr, "Missing length of read by default: 100\n");
			taille_fenetre = 100;
		}
		else
		{
			if (taille_fenetre==-1)
			{
				// fprintf(stderr, "length of the read = length sequence\n");
				taille_fenetre = taille_seq;
			}
		}
		
		if (taille_fenetre > taille_seq)
		{
			fprintf(stderr, "the read size is greater than that of the sequence %d > %d (seq n°%d)\n",taille_fenetre,taille_seq,x);
		}


			nombre_sous_sequence = taille_seq - taille_fenetre+1;
					// printf("nombre_sous_sequence = %d\n",nombre_sous_sequence);

		// printf("Nombre seq = %d\n",nombre_sous_sequence);
		/******** ALLOC ************/

		// char **accession_tab = (char **)(malloc(MAX_SEQ*sizeof(char *)));
		// my_error (accession_tab,"Erreur malloc accession_tab");

		// for (i=0; i < MAX_SEQ; i++)
		// {
		// 	accession_tab[i] = (char *)(malloc(20*sizeof(char)));
		// 	my_error(accession_tab[i],"Erreur malloc accession_tab[i]");
		// }
		// resultat = (int ***)(malloc(nb_kmer*sizeof(int **)));
		
		resultat = (int ***)(calloc(nb_kmer,sizeof(int **)));
		my_error(resultat,"Erreur malloc resultat");
	
		
		

		for (i = 0; i < nb_kmer; i++)
		{
			// resultat[i] = (int **)(malloc(nombre_sous_sequence*sizeof(int *)));
			
			resultat[i] = (int **)(calloc(nombre_sous_sequence,sizeof(int *)));
			my_error(resultat,"Erreur malloc resultat[i]");
			

			nb_possibilite = puissance4(taille_kmer[i]);
			// printf("nb_possibilite = %d\n",nb_possibilite);


			for(j=0; j < nombre_sous_sequence;j++)
			{
				// printf("Debut Alloc interne\n");	
				// printf("i = %d \t j = %d\n",i,j);	
				// resultat[i][j] = (int *)(malloc(nb_possibilite*sizeof(int )));
				resultat[i][j] = (int *)(calloc(nb_possibilite,sizeof(int)));
				// printf("Fin Alloc interne\n");
				my_error(resultat,"Erreur malloc resultat[i][j]");
				// accession_tab[indice_ligne++] = les_accessions[x];
			}

		}
		/********** ALLOC FIN **********/


		// int *tab_comptage = calloc((taille+ 1), sizeof(int));


		/**********************************
		**
		** Routine principale
		***********************************/

		for(i=0;i<nb_kmer;i++)
		{
			if (taille_fenetre < taille_pattern[i])
			{	
				// printf("taille_fenetre_origine = %d et j = %d\n", taille_fenetre_origine,j);
				fprintf(stderr, "the pattern size is greater than that of the read\n");
				exit(EXIT_FAILURE);
			}


			count(les_sequences[x],taille_seq,tab_bool[i],&taille_kmer[i],resultat[i],taille_fenetre);

		}



		// imprime_csv(out,nb_kmer,taille_kmer,nombre_sous_sequence,resultat);

		imprime_weka(out,nb_kmer,taille_kmer,nombre_sous_sequence,resultat,taxid,x,les_accessions[x]);

	}


	/**************************
	**********FREE*************
	***************************/
	free(taille_kmer);
	free(taille_pattern);
	// free(tab_bool);
	free_tab_bool(tab_bool,15);

	
	// printf("Debut Free\n");
	for(i=0; i < nb_kmer; i++)
	{
		for(j=0; j < nombre_sous_sequence; j++)
		{
			// printf("Debut free interne\n");	
			// printf("i = %d \t j = %d\n",i,j);	
			free(resultat[i][j]);
			// printf("Fin free interne\n");	
		}
		// printf("DEBUT FREE\n");
		
		// printf("FIN FREE\n");
	}
	// printf("FIN FREE\n");
	for(i=0; i < nb_kmer; i++)
	{
		
		free(resultat[i]);
		
	}


	free(resultat);
	
	les_sequences[1] = (char *)(calloc(TAILLE_MAX_SEQ,sizeof(char)));
	les_sequences[2] = (char *)(calloc(TAILLE_MAX_SEQ,sizeof(char)));

	for (i = 0; i < MAX_SEQ; i++)
	{
		printf("Je vais liberer les_sequences[%d]\n",i);
		// if (i==1 || i==0)
			free(les_sequences[i]);
		printf("Fin liberation\n");
		// free(les_accessions[i]);
	}

	for (i = 0; i < MAX_SEQ; i++)
	{
		printf("Je vais liberer les_accessions[%d]\n",i);
		if (i==1 || i==0)
			free(les_accessions[i]);
		printf("Fin liberation\n");
		// free(les_accessions[i]);
	}


	free(les_sequences);
	free(les_accessions);



	return 0;
	
}
