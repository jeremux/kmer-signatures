#!/usr/bin/perl -w


use strict;

use Getopt::Long;
use Cwd 'abs_path';

my $help;
my $fichier_pattern;
my $dossier;
my $nom_sequence;
my $fichier_feuille;
my $taille_read = -1;
my %hash_taxid = ();

my @taille_kmer;
my $nb_kmer = 0;



GetOptions ('pattern|p=s' => \$fichier_pattern,
            'read|r=i' => \$taille_read,	
            'leaf|f=s' => \$fichier_feuille,	
            'help|h' => \$help);


print "fichier_feuille = $fichier_feuille\n";
open (FEUILLE,'<',$fichier_feuille) || die "Can't open file $fichier_feuille:$!\n";
open (KMER,  '<', $fichier_pattern) || die "Can't open file $fichier_pattern:$!\n";

my $premiere_ligne = <FEUILLE>;

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
 

# print WEKA "\@relation kmots\n\n";
my $tmp_size = 0;

# foreach my $k (@taille_kmer)
# {
# 	&write_entete($k);
# }

# print WEKA "\@attribute ID {";

my $cpt;

# foreach my $k (@les_taxids)
# {
# 	if(++$cpt == scalar(@les_taxids))
# 	{
# 		print WEKA "$k}\n\n\@data\n";
# 	}
# 	else
# 	{
# 		print WEKA "$k,";
# 	}
# }

# ROUTINE COMPTAGE

while (<FEUILLE>)
{
	my $k = $_ ;
	$k =~ s/^\s+//;
	$k =~ s/\s+$//;
	# print "k = $k\n";

	system("mkdir -p $k/count");
	my $les_genomes = `ls $k/genomes_*.fasta`;

	my @liste_genomes = split('\n',$les_genomes);

	my $nom_fichier = $k.'/count/count_S'.$taille_read.'.arff';


	open(WEKA,'>',$k.'/count/count_S'.$taille_read.'.arff') || die "Can't open file $nom_fichier:$!\n";
	# print "*****les_genomes*****\n";
	foreach my $j (@liste_genomes)
	{
		# print "nom_fichier = $nom_fichier\n";
		# system("./../count_kmer/src/count_kmer -i $j -k $fichier_pattern -l $taille_read -o $nom_fichier -t X");
		# print "./../count_kmer/src/count_kmer -i $j -k $fichier_pattern -l $taille_read -o $nom_fichier -t X >> $nom_fichier\n";
		system("./../count_kmer/src/count_kmer -i $j -k $fichier_pattern -l $taille_read -o $nom_fichier -t X >> $nom_fichier");

	}
	close(WEKA);
	# 
}

# system("rm A954xRwm");
# print WEKA "@relation "





