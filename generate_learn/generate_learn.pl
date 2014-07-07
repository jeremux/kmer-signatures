#!/usr/bin/perl -w


use strict;

use Getopt::Long;
use Cwd 'abs_path';
use Parallel::ForkManager;
use threads;
use threads::shared;

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


#on ne génère pas d'apprentissage
#si on est a une feuille
#en l'occurence si ls -d est un echec
my $sous_dossier = `ls -d $dossier/*/ 2> A954xRwm`;
	# print "valeur = $?\n";
if ($? != 0) {
	exit;
}

#les taxids des dossiers
#du sous dossier courant
my @les_taxids;

print "sous_dossier = $sous_dossier\n";
#les sous dossiers
#du dossier courant
my @les_sous_dossier = split('\n',$sous_dossier);

#tableau des paths des 
# sous dossiers
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

#en tete weka avec une taille 
#de k-mer
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

	#impression des k-mer pour weka
	foreach my $w (@words)
	{
		print WEKA "\@attribute $w numeric\n";
	}
}
 

my $nom_fichier = $path_racine.'/learn/learning_S'.$taille_read.'.arff';


open(WEKA,'>',$path_racine.'/learn/learning_S'.$taille_read.'.arff') || die "Can't open file $nom_fichier:$!\n";

###
#premiere ligne weka
print WEKA "\@relation kmots\n\n";
my $tmp_size = 0;

#Pour chaque k-mer
#generer la liste des k-mers
foreach my $k (@taille_kmer)
{
	&write_entete($k);
}

#les attributs à prédire
print WEKA "\@attribute ID {";

my $cpt;

#pour chaque taxid
foreach my $k (@les_taxids)
{
	#si c'est le dernier de la liste
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

#pour chaque sous dossier
#on effectue le comptage

my $nb_sous_dossier = @path_sous_dossier;
my @childs;
my @threads;
my $count = 0;

sub doThread
{
	for ( my $count = 1; $count <= $nb_sous_dossier; $count++) {
			my $tmp = $count - 1;
			my $k = $path_sous_dossier[$tmp];
	        my $t = threads->new(\&process, $count,$k);
	        push(@threads,$t);
	}
	foreach (@threads) {
	        my $num = $_->join;
	        print "done with $num\n";
	}
}

sub doFork
{
	for ( my $count = 1; $count <= $nb_sous_dossier; $count++) 
	{
	        my $pid = fork();
	        my $tmp = $count - 1;
	        my $k = $path_sous_dossier[$tmp];
	        if ($pid) {
	        # parent
	        # print "pid is $pid, parent $$\n";
	        push(@childs, $pid);
	        } elsif ($pid == 0) {
	                # child
	                &process($count,$k);
	                exit 0;
	        } else {
	                die "couldnt fork: $!\n";
	        }
	}


	foreach (@childs) {
	        my $tmp = waitpid($_, 0);
	         # print "done with pid $tmp\n";
	}
}

sub do_nothing
{
	for ( my $count = 1; $count <= $nb_sous_dossier; $count++) 
	{
	        my $tmp = $count - 1;
	        my $k = $path_sous_dossier[$tmp];
	       &process($tmp,$k);
   }
}

sub process
{
	my ($num,$k) = (@_);
	my $txid = $hash_taxid{$k};
	print "txid = $txid\n";
	my $les_genomes = `ls $k/genomes_*.fasta`;
	my @liste_genomes = split('\n',$les_genomes);
	print "*****les_genomes*****\n";
	foreach my $j (@liste_genomes)
	{
		print "******COMMANDE*********\n";
		print "./../count_kmer/src/count_kmer -i $j -k $fichier_pattern -l $taille_read -o $nom_fichier -t $txid\n";
		print "************************\n";
		system("./../count_kmer/src/count_kmer -i $j -k $fichier_pattern -l $taille_read -o $nom_fichier -t $txid");

	}
	sleep $num;
    # print "done with child process for $num\n";
    return $num;
	# 
}

# &doFork();
# &doThread();
&do_nothing();

close(WEKA);

#suppression du fichier temporaire d'erreur
system("rm A954xRwm");
# print WEKA "@relation "





