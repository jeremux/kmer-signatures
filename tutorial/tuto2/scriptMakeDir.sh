#!/bin/sh
cd ../create_db/;

#########################################################################
./generateDirectories.pl -id 2759 -gen ../qwery/eukaryota.gb -bound 10 ##
#########################################################################

if [ $? -ne 0 ]
then
   echo "error "
   exit 3 
fi

echo "../create_db/Eukaryota__2759 generated"
echo "../generate_data/generateGenbank_Eukaryota.sh generated"
echo "../generate_data/listGenbank.txt generated"
echo "Now you can go to the tuto3 !"
