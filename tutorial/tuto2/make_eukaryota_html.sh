#!/bin/sh

cd ../create_db;
ktImportXML krona_Eukaryota.xml -o eukaryota.html;

echo "../create_db/eukaryota.html generated";
