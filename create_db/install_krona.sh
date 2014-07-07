#!/bin/sh

echo ""
echo ""
echo "*********************"
echo "Debut install krona"
echo "*********************"
echo ""
echo ""

if [ -d "krona" ]; then
	echo "Le dossier krona existe"
	exit 1;
else

git clone git://git.code.sf.net/p/krona/git krona;
cd krona/KronaTools
sudo perl install.pl
cd ../..
# rm -rf ./n66d5Ey8pBBHb6

echo ""
echo ""
echo "**************************"
echo "Fin install krona, exemple"
echo "**************************"
echo ""
echo ""
echo "*********************"
echo "ktImportXML krona_Eukaryota.xml"
echo "*********************"
echo ""
echo ""
fi


