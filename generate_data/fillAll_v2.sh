#!/bin/bash

cp feuille.txt feuille2.txt;

racine=`head -n 1 feuille2.txt`;
racine="${racine%?}";

# suppression de la racine
sed -i -e "1d" feuille2.txt;


cpt=0;

#pour chaque dossier feuille
for line in $(cat feuille2.txt); do
	cpt=$((cpt+1));
	#on bouge dans le dossier
	cd $line;
	#pour chaque chemin des fichiers (-tyep f) de ce dossier
	for file in `find \`pwd\` -maxdepth 1 -type f`; do	
		#on recupere le nom du fichier
		nom_fichier=$(basename "$file");
		
		#tant qu'on est pas a la racine, on fait lien seymbolique
		cd ..;
		courant=`pwd`;
		while [ "$courant" != "$racine" ] 
		do
			# ln -s $file $nom_destination;
			ln -s $file $nom_fichier;
			cd ..;
			courant=`pwd`
		done
		#on rebouge dans le dossier
		cd $line;
	done
done

rm -rf feuille2.txt;