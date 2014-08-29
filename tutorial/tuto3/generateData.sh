#!/bin/sh

cd ../generate_data;

#################################
./generateGenbank_Eukaryota.sh; #
#################################

if [ $? -ne 0 ]
then
   echo "error "
   exit 3 
fi

echo "All data in ../create_db/Eukaryota__2759 generated"
echo "Ready for the 4th tuto ?"
