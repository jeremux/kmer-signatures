#!/bin/sh
cd ../create_db;
sh install_krona.sh;

if [ $? -ne 0 ]
then
   echo "error "
   exit 3 
fi

echo "install krona in ../create_db/ done"
