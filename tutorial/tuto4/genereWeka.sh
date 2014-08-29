#!/bin/sh

cd ../cpp/count_kmer
make realclean
make

############################################################################
./count_kmer -f seqTuto.fasta -k patternTuto.txt -o tuto.arff -l -1 --weka #
############################################################################


if [ $? -ne 0 ]
then
   echo "error "
   exit 3 
fi

make realclean

echo "../cpp/count_kmer/tuto.arff generated"
