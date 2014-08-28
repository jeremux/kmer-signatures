#!/bin/bash
echo "Traitement 3/3 (linkage)";
origine=`pwd`;
cp feuille.txt feuille2.txt;

racine=`head -n 1 feuille2.txt`;
racine="${racine%?}";

# suppression de la racine
sed -i -e "1d" feuille2.txt;

pathAAcox1="data/fasta/aminoAcids/cox1/";
pathAAcox2="data/fasta/aminoAcids/cox3/";
pathAAcox3="data/fasta/aminoAcids/cox2/";
pathAAcytb="data/fasta/aminoAcids/cytb/";

pathNTcox1="data/fasta/nucleotides/cox1/";
pathNTcox2="data/fasta/nucleotides/cox2/";
pathNTcox3="data/fasta/nucleotides/cox3/";
pathNTcytb="data/fasta/nucleotides/cytb/";
pathNTgen="data/fasta/nucleotides/genomes/";

pathGenbank="data/genbank/";


cpt=0;

array=("$pathNTcox1" "$pathNTcox2" "$pathNTcox3" "$pathNTcytb" 
	   "$pathAAcox1" "$pathAAcox2" "$pathAAcox3" "$pathAAcytb"
	   "$pathGenbank" "$pathNTgen"
	  )



#routineLN "Ciliophora__5878" "data/fasta/aminoAcids/cox1"
routineLN()
{
	for file in `find \`pwd\` -maxdepth 1 -type f`; do	
		#on recupere le nom du fichier
		nom_fichier=$(basename "$file");
		#on rebouge doans lespece
		cd $1;
		#tant qu'on est pas a la racine, on fait lien symbolique
		cd ..;
		courant=`pwd`;
		while [ "$courant" != "$racine" ] 
		do
			# ln -s $file $nom_destination;
			tmp=`pwd`;
			mkdir -p $2;
			ln -s $file $2/$nom_fichier;
			cd ..;
			courant=`pwd`
		done
		#on rebouge dans le dossier
		cd $line;
	done
}


routineMain()
{
	# pour chaque dossier feuille
	for line in $(cat feuille2.txt); do
		cpt=$((cpt+1));
		#on bouge dans le dossier
		cd $line;
		for i in "${array[@]}"; do 
			# lieu=`pwd`;
			# echo "je suis dans $lieu"; 
    		cd $i;
    		routineLN $line $i;
		done
		#pour chaque chemin des fichiers (-tyep f) de ce dossier

	done
}

routineMain;

cd $origine;
rm -rf feuille2.txt;