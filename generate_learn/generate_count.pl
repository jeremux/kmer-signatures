#!/usr/bin/perl -w


use strict;
use Parallel::ForkManager;
use threads;
use threads::shared;
use Getopt::Long;
use Cwd 'abs_path';

my $help;
my $fichier_pattern;
my $dossier;
my $nom_sequence;
my $fichier_feuille;
my $taille_read = -1;
my %hash_taxid = ();
my $nb_thread = 0;
my $p = 0;

my @taille_kmer;
my $nb_kmer = 0;

# En option le fichier contenant la graine, la taille du read, 
# et les feuilles (liste de path) de la bdd. si p = nb de fork
GetOptions ('pattern|p=s' => \$fichier_pattern,
            'read|r=i' => \$taille_read,	
            'leaf|f=s' => \$fichier_feuille,
            'thread|t=i' => \$p,	
            'help|h' => \$help);

# on vérifie si le programme count_kmer a bien
# été compilé
if (! -e "./../count_kmer/src/count_kmer") 
{
	print "Lancez le makefile du dossier ./../count_kmer/src\n";
	exit 1;
} 

# print "fichier_feuille = $fichier_feuille\n";
open (FEUILLE,'<',$fichier_feuille) || die "Can't open file $fichier_feuille:$!\n";
open (KMER,  '<', $fichier_pattern) || die "Can't open file $fichier_pattern:$!\n";

# la premiere ligne du fichier feuille est le path de la bdd
my $premiere_ligne = <FEUILLE>;

my $tmp_size = 0;

my $cpt;

# ROUTINE COMPTAGE

my %hash_chemins = ();


my $nb_feuilles = 0;
while (<FEUILLE>)
{
	my $k = $_ ;
	#trim des whitespaces
	$k =~ s/^\s+//;
	$k =~ s/\s+$//;

	$hash_chemins{1} .= "," . $k ;
	$nb_feuilles++;

	#split des path en plusieurs parties
	for (my $i = 2 ; $i <= $nb_thread ;$i++)
	{
		$k = <FEUILLE> ;
	
	
		if (defined($k))
		{
			$k =~ s/^\s+//;
			$k =~ s/\s+$//;
			# print "Je vais push $k\n";
			$hash_chemins{$i} .= "," . $k ;
			$nb_feuilles++;
		}
	}
}

my @childs;
my @threads;
my $count = 0;


########################################################
########################################################
###############    FORK ################################
########################################################
########################################################

sub doThread
{
	for ( my $count = 1; $count <= $nb_thread ; $count++) 
	{
	        my $pid = fork();
	        my $chemin = $hash_chemins{$count};
	        if ($pid) {
	        # parent
	        # print "pid is $pid, parent $$\n";
	        push(@childs, $pid);
	        } elsif ($pid == 0) {
	                # child
	                &writeThread($count,$chemin);
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

sub normal
{
	my ($h) = (@_);
	my $chemin = $hash_chemins{1};
	&writeNormal($chemin);
}

sub main
{
	if ($p==1)
	{
		&doThread();
	}
	else
	{
		&normal(%hash_chemins);
	}
}

sub count_routine
{
	my ($chemins) = (@_);
	my @tab_chemins = split(',',$chemins);

        foreach my $k (@tab_chemins)
        {
        	if ($k ne "")
        	{
		        # print "started child process for  $num\n";
		        system("mkdir -p $k/count");
				my $les_genomes = `ls $k/genomes_*.fasta`;

				my @liste_genomes = split('\n',$les_genomes);

				my $nom_fichier = $k.'/count/count_S'.$taille_read.'.arff';


				open(WEKA,'>',$k.'/count/count_S'.$taille_read.'.arff') || die "Can't open file $nom_fichier:$!\n";
				# print "*****les_genomes*****\n";
				foreach my $j (@liste_genomes)
				{
					print "###########################\n";
					print "./../count_kmer/src/count_kmer -i $j -k $fichier_pattern -l $taille_read -o $nom_fichier -t X\n";
					system("./../count_kmer/src/count_kmer -i $j -k $fichier_pattern -l $taille_read -o $nom_fichier -t X ");
					#system("./../count_kmer/src/count_kmer -i $j -k $fichier_pattern -l $taille_read -o $nom_fichier -t X >> $nom_fichier");
					print "###########END#############\n\n";

				}
				close(WEKA);
			}
		}
}
 
sub writeThread {
        my ($num,$chemins) = (@_);

        # print "CHEMIN = $chemins \n";
        &count_routine($chemins);
        sleep $num;
        # print "done with child process for $num\n";
        return $num;
}

sub writeNormal {
        my ($chemins) = (@_);

        my @tab_chemins = split(',',$chemins);

        &count_routine($chemins);
}


&main();



# system("rm A954xRwm");
# print WEKA "@relation "





