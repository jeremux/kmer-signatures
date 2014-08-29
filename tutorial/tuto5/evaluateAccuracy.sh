#!/bin/sh

cd ../cpp/count_kmer

./compilBayesJava.sh

echo "=========="
echo "Learn size = -1 VS predict size = 100" 

echo -n "Accuracy = "

###############################################################################################
sh execCrossVal.sh "../../create_db/Eukaryota__2759/Alveolata__33630/" "10" "100" "complete" ##
###############################################################################################

if [ $? -ne 0 ]
then
   echo "error "
   exit 3 
fi

echo ""

echo "Now you are ready for the tuto 6 :)"
