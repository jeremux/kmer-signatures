#!/usr/bin/perl -w


use strict;

use Getopt::Long;
use Cwd 'abs_path';

my $help;
my $fichier_pattern;
my $dossier;
my $nom_sequence;
my $taille_read = -1;
my %hash_taxid = ();

my @taille_kmer;
my $nb_kmer = 0;



GetOptions ('dir|d=s'	 => \$dossier,
            'pattern|p=s' => \$fichier_pattern,
            'read|r=i' => \$taille_read,		
            'help|h' => \$help);


my $sous_dossier;

eval {
	my $sous_dossier = `ls -d $dossier/*/ 2> A954xRwm`;
	# print "valeur = $?\n";
};
if ($? != 0) {
	exit;
}

my @les_taxids;

# print "total = $sous_dossier";

my @les_sous_dossier = split('\n',$sous_dossier);
my @path_sous_dossier;
my $others = "others";
my $path_racine = abs_path($dossier);

print "path_racine = $path_racine\n";
system("mkdir -p $path_racine/learn/");


# je récupère mes taxid
foreach my $d (@les_sous_dossier)
{
	# print "d = $d\n";
	if ($d =~ m/(.*)__(.*)\// )
	{
		my $tmp = abs_path($d);
		my $abs = abs_path($d);
		if (!($d =~ /learn/))
		{
			push(@path_sous_dossier,$abs);
			if ($d =~ /$others/)
			{
				$hash_taxid{$abs} = "others";
				push(@les_taxids, "others");
			}
			else
			{
				$hash_taxid{$abs} = "$2";
				push(@les_taxids, "$2");
			}	
		}
		
	}
}

foreach my $x (@path_sous_dossier) 
{
	print "path_sous_dossier = $x\n";
}

open (KMER,  '<', $fichier_pattern) || die "Can't open file $fichier_pattern:$!\n";

while (<KMER>) 
{
	my $l = $_;
	#trim
	$l =~ s/\d//g;
	$l =~ s/\s*//g;	

	my @tab_kmer = split(//,$l);
	my $tmp_size = 0;

	foreach my $k (@tab_kmer)
	{
		if($k eq "#")
		{
			$tmp_size++;
		}

	}

	if ($tmp_size != 0)
	{
		print "taille_kmer = $tmp_size\n";
		push(@taille_kmer,$tmp_size);
		$nb_kmer++;
	}
	
}

sub write_entete
{

	my ($k) = @_;
	my @bases = ('A','C','G','T');
	my @words = @bases;
	my @newwords;
	for my $i (1..$k-1)
	{
		undef @newwords;
		foreach my $w (@words)
		{
			foreach my $b (@bases)
			{
				push (@newwords,$w.$b);
			}
		}
		undef @words;
		@words = @newwords;
	}
	foreach my $w (@words)
	{
		print WEKA "\@attribute $w numeric\n";
	}
}
 

my $nom_fichier = $path_racine.'/learn/learning_S'.$taille_read.'.arff';


open(WEKA,'>',$path_racine.'/learn/learning_S'.$taille_read.'.arff') || die "Can't open file $nom_fichier:$!\n";

print WEKA "\@relation kmots\n\n";
my $tmp_size = 0;

foreach my $k (@taille_kmer)
{
	&write_entete($k);
}

print WEKA "\@attribute ID {";

my $cpt;

foreach my $k (@les_taxids)
{
	if(++$cpt == scalar(@les_taxids))
	{
		print WEKA "$k}\n\n\@data\n";
	}
	else
	{
		print WEKA "$k,";
	}
}

# ROUTINE COMPTAGE

foreach my $k (@path_sous_dossier)
{
	my $txid = $hash_taxid{$k};
	print "txid = $txid\n";
	my $les_genomes = `ls $k/genomes_*.fasta`;
	my @liste_genomes = split('\n',$les_genomes);
	print "*****les_genomes*****\n";
	foreach my $j (@liste_genomes)
	{
		# system("./../count_kmer/src/count_kmer -i $j -k $fichier_pattern -l $taille_read -o $nom_fichier -t $txid");
		system("./../count_kmer/src/count_kmer -i $j -k $fichier_pattern -l $taille_read -o $nom_fichier -t $txid >> $nom_fichier");

	}

	# 
}

close(WEKA);
system("rm A954xRwm");
# print WEKA "@relation "





