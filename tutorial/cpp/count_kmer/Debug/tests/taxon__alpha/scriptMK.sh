#!/bin/bash

pathAAcox1="data/fasta/aminoAcids/cox1/";
pathAAcox2="data/fasta/aminoAcids/cox3/";
pathAAcox3="data/fasta/aminoAcids/cox2/";
pathAAcytb="data/fasta/aminoAcids/cytb/";

pathNTcox1="data/fasta/nucleotides/cox1/";
pathNTcox2="data/fasta/nucleotides/cox2/";
pathNTcox3="data/fasta/nucleotides/cox3/";
pathNTcytb="data/fasta/nucleotides/cytb/";
pathNTgen="data/fasta/nucleotides/genomes/";

pathGenbank="data/genbank/";

array=("$pathNTcox1" "$pathNTcox2" "$pathNTcox3" "$pathNTcytb" 
	   "$pathAAcox1" "$pathAAcox2" "$pathAAcox3" "$pathAAcytb"
	   "$pathGenbank" "$pathNTgen"
	  )

for i in "${array[@]}"; do 
	mkdir -p $i;
done


