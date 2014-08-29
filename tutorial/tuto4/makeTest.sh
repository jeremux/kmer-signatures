#!/bin/sh

cd ../cpp/count_kmer
make realclean
make
##################
./count_kmer -t ##
##################
if [ $? -ne 0 ]
then
   echo "error "
   exit 3 
fi

##################
./count_kmer -T ##
##################
if [ $? -ne 0 ]
then
   echo "error "
   exit 3 
fi

echo "test done"
