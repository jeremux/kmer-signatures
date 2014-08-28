#!/bin/sh

cd ../generate_graph

perl make_graph.pl -in ../cpp/count_kmer/result.log -title tutoCox1

echo "../generate_graph/tutoCox1.pdf generated"

echo "Now tuto8 ? No it's the end, enjoy !"
