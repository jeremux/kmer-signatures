#!/bin/bash

rm -rf leaf.txt;

recursiveSearch() {
    [ "`ls "$1" | wc -l`" == "0" ] && echo "$1" >> leaf.txt # Si le dossier est vide
    for file in "$1"/*; do
        if [ -d "$file" ]; then
            recursiveSearch "$file"
        fi
    done
}


# echo "recherche dans $1"
recursiveSearch $1;

#on enleve le double slash du debut des paths
sed -i -e "s/\/\//\//g" leaf.txt;

rm -rf ../generate_data/listGenbank.txt ;
rm -rf scrpit_clean_$3.sh ;

echo "find $1 -type f -exec rm '{}' \;" >> scrpit_clean_$3.sh;
echo "find $1 -type l -exec rm '{}' \;" >> scrpit_clean_$3.sh;

for line in $(cat leaf.txt); do 
	cat $2 | grep $line >> ../generate_data/listGenbank.txt ;
done

rm -rf leaf.txt
chmod +x scrpit_clean_$3.sh;