#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>



const unsigned char m[256] =
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

int **make_tab_bool(int nb_kmer,int taille_max_kmer)
{
	int i;

	int **les_kmers = (int **)calloc(nb_kmer,(sizeof(int *)));

	if(les_kmers==NULL)
	{
		fprintf(stderr, "Erreur creation tableau\n");
		exit(EXIT_FAILURE);
	}

	for(i=0;i<nb_kmer;i++)
	{
		les_kmers[i] = (int *)calloc(taille_max_kmer,(sizeof(int)));
		if (les_kmers[i]==NULL)
		{
			fprintf(stderr, "Erreur creation tableau\n");
			exit(EXIT_FAILURE);
		}
	}

	return les_kmers;
}

void free_tab_bool(int **tab,int nb_kmer)
{
	int i;
	for(i=0;i<nb_kmer;i++)
	{
		free(tab[i]);
	}

	free(tab);
}

void free_init(int m, char **sequences)
{
	int i;
	for (i = 0; i < m; i++)
	{
		sequences[i] = (char *)(calloc(m,sizeof(char)));
	}
}

int read_pattern(const char *path,int **tab_bool)
{
	int c,val;
	int c_last;
	int indice_kmer = 0;
	int indice_carac = 0;

	FILE* fichier = NULL;

	fichier = fopen(path,"r");

	if (fichier != NULL)
	{
		while((c = fgetc(fichier)) != EOF)
		{
			// printf("Valeur de c = %d\n",c);
			if (c=='\n')
			{
				if (c!=c_last)
				{
					// printf("indice_carac = %d\n",indice_carac);
					tab_bool[indice_kmer][indice_carac] = -1;
					indice_kmer++;
				}
				indice_carac = 0;
				
			}	
			else
			{
				if (c!=' ')
				{
					if (c=='#')
					{
						val = 1;
					}
					if (c=='_')
					{
						val = 0;
					}

					tab_bool[indice_kmer][indice_carac] = val;
					indice_carac++;
				}
			}

			c_last = c;
		}

		if(c_last=='#' || c_last=='_')
		{
			tab_bool[indice_kmer][indice_carac] = -1;
			indice_kmer++;
		}

		fclose(fichier);
	}
	else
	{
		fprintf(stderr,"Erreur ouverture %s\n",path);
	}

	// printf("tab_bool3[%d][%d] = %d\n",0,0,tab_bool[0][0]);
	// printf("Nombre de kmer = %d\n",indice_kmer);
	// printf("Taille du kmer = %d\n",indice_carac);
	// printf("tab_bool[%d] = %d\n",indice_carac,tab_bool[0][indice_carac+1]);


	return indice_kmer;
}

void append(char* s, char c)
{
        int len = strlen(s);
        s[len] = c;
        s[len+1] = '\0';
}

char **read_seq(const char *path,int *nb_sequences,char **les_sequences, char **les_accessions)
{

	int c,i;
	int flag_titre = 0;
	int flag_take_acc = 0;
	int flag_take_seq = 0;
	int ancien_char1='\0';
	int ancien_char2='\0';
	int x = -1;
	int flag_tmp = 0;
	FILE* fichier = NULL;
	int max = 0;

	// char *result = (char *)(calloc(TAILLE_MAX_SEQ,sizeof(char)));
	// result[0] = '\0';

	i = -1;
	

	fichier = fopen(path,"r");

	if (fichier != NULL)
	{
		while((c = fgetc(fichier)) != EOF)
		{
			// printf("carac = %c\n",c);
			if (c=='>')
			{
				flag_titre = 1;
				flag_take_seq = 0;
				goto end;
			}

			if (flag_titre && ancien_char2=='_' && ancien_char1=='_')
			{	
				flag_take_acc = 1;
				i++;
				if (i>0)
				{
					les_sequences[i-1][x]='\0';
				}
				x = 0;
			}

			// printf("test\n");

			if (c=='\n' && flag_titre)
			{
				flag_take_seq = 1;
				flag_titre = 0;
				flag_take_acc = 0;
				// printf("\n");
				// i++;

				// if(i>=0)
				// {
				// 	printf("Im gonna print\n");
				// 	printf("les_accessions[%d] = %s\n",i,les_accessions[i]);
				// 	if( i>0)
				// 	{
				// 		printf("\n******************\n");
				// 		printf("les_sequences[%d] = %s\n",i-1,les_sequences[i-1]);
				// 	}
				// }
				// printf("Lecture les_sequence %d\n",i);
			}

			if(flag_take_acc)
			{
				// printf("%c",c);
				// printf("Ajout de %c à les_accessions[%d]\n",c,i);
				append(les_accessions[i],c);
			}

			if (flag_take_seq && c!=' ' && c!='\n' && c!='\t')
			{
				flag_tmp = 0;
				switch(c)
				{

					case 'y':
						c = 'c';
						break;
					case 's':
						c = 'g';
						break;
					case 'w':
						printf("rencontre\n");
						flag_tmp = 1;
						c = 't';
						break;
					case 'k':
						c = 'a';
						break;
					case 'm':
						c = 'c';
						break;
					case 'b':
						c = 'g';
						break;
					case 'd':
						c = 'a';
						break;
					case 'h':
						c = 't';
						break;
					case 'v':
						c = 'c';
						break;
					case 'n':
						c = 'a';
						break;
					default:
						break;
				}
				// append(les_sequences[i],c);
				if(x >= max)
				{
					max = x;
				}
				
				les_sequences[i][x++]=c;
				if (flag_tmp) printf("ajout de %c dans %d\n",c,i);
			}

			ancien_char2 = ancien_char1;
			ancien_char1 = c;

			end: ;
		}

		// printf("%s\n",les_sequences[i]);
		les_sequences[i][x]='\0';
		*nb_sequences = i + 1;
		// printf("max = %d\n",max);
		fclose(fichier);

		// printf("\n******************\n");
		// printf("les_sequences[%d] = %s\n",i,les_sequences[i]);
	}
	else
	{
		fprintf(stderr,"Erreur ouverture %s\n",path);
	}

	// printf("tab_bool3[%d][%d] = %d\n",0,0,tab_bool[0][0]);
	// printf("Nombre de kmer = %d\n",indice_kmer);
	// printf("Taille du kmer = %d\n",indice_carac);
	// printf("tab_bool[%d] = %d\n",indice_carac,tab_bool[0][indice_carac+1]);
	// printf("final %s\n",result);
	return les_sequences;
}