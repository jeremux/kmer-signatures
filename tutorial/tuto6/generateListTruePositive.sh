#!/bin/sh

cd ../cpp/count_kmer

perl execAndEval.pl --root "../../create_db/Eukaryota__2759/Alveolata__33630" --start "100" --end "350" --step "50" --list list.txt --kmer "patternTuto.txt" --sample "20" --key "cox1"

echo "../cpp/count_kmer/result.log generated"

echo "================="
echo "Now, the last tutorial: tuto7"
