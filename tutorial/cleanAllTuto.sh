#!/bin/sh

for f in tuto1 tuto2 tuto3 tuto4 tuto5 tuto6 tuto7 ; do

	cd $f;
	
	sh scriptClean.sh 2> /dev/null
	echo "*********************"
	echo "*** $f cleaned ***"
	echo "*********************"
	cd ..
done
