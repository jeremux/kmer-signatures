#!/bin/sh

cd ../cpp/count_kmer
make realclean
make
./count_kmer -t
./count_kmer -T

echo "test done"
