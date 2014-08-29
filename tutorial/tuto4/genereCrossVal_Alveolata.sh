#!/bin/sh

cd ../cpp/count_kmer
make realclean
make
################################################################################################################################################################################
./count_kmer --sample "20" --start "100" --step "50" --end "350" --learn "-1" --kmer "patternTuto.txt" --root "../../create_db/Eukaryota__2759/Alveolata__33630/" --key "cox1" #
################################################################################################################################################################################
if [ $? -ne 0 ]
then
   echo "error "
   exit 3 
fi


make realclean

echo "weka data generated in ../../create_db/Eukaryota__2759/Alveolata__33630/frequencies/"

echo "==========="
echo "Next step ?"
