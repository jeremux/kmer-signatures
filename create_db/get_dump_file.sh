#!/bin/sh

wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz ;
mkdir -p bdd ;
tar -C bdd/ -xzf taxdump.tar.gz ;
rm -rf taxdump.tar.gz ;