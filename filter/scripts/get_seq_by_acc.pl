#!/usr/bin/perl -w


use strict;
use warnings;

# Chargement du module GenBank
use Bio::DB::GenBank;

# Création du handle permettant de se connecter à la banque de données GenBank avec le constructeur new qui ne prend aucun argument
my $gb = new Bio::DB::GenBank;
my $file = shift @ARGV;
my $seq1 = "";
my $description = "";
my $acc = "";

# Création de l'objet $seq1 correspondant à l'accession J00522
open (FICHIER_IN, '<', $file) || die "Can't open file:$!\n"; 

while(<FICHIER_IN>)
{
	$acc = $_;
	chomp($acc);

	if ($acc ne "")
	{
		$seq1 = $gb->get_Seq_by_acc($acc);
		$description = $seq1->desc();
	}

	print "$description\n";
}

