#!/bin/sh

echo ""
echo ""
echo "*********************"
echo "Debut install krona"
echo "*********************"
echo ""
echo ""

if [ -d "n66d5Ey8pBBHb6" ]; then
	echo "Le dossier n66d5Ey8pBBHb6 existe"
	exit 1;
else

git clone git://git.code.sf.net/p/krona/git n66d5Ey8pBBHb6;
cd n66d5Ey8pBBHb6/KronaTools
sudo perl install.pl
cd ../..
rm -rf ./n66d5Ey8pBBHb6

echo ""
echo ""
echo "**************************"
echo "Fin install krona, exemple"
echo "**************************"
echo ""
echo ""
echo "*********************"
echo "ktImportXML fichier_krona.xml"
echo "*********************"
echo ""
echo ""
fi


