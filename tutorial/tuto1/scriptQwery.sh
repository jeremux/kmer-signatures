#!/bin/sh
cd ../qwery ;

#############################################################################
perl qweryNCBI.pl -id 2759 -m toto@mail.com -out eukaryota -more "NOT HOMO" #
#############################################################################

if [ $? -ne 0 ]
then
   echo "error "
   exit 3 
fi

echo "../qwery/eukaryota.gb downloaded"
echo "Now you can go to the tuto2 !"
