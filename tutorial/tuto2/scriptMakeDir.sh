#!/bin/sh
cd ../create_db/;
./generateDirectories.pl -id 2759 -gen ../qwery/eukaryota.gb -bound 10
echo "../create_db/Eukaryota__2759 generated"
echo "../generate_data/generateGenbank_Eukaryota.sh generated"
echo "../generate_data/listGenbank.txt generated"
echo "Now you can go to the tuto3 !"
