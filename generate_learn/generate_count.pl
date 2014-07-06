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

my %hash_chemins = ();
my @les_chemins1;
my @les_chemins2;
my @les_chemins3;
my @les_chemins4;



my $nb_feuilles = 0;
while (<FEUILLE>)
{
	my $k = $_ ;
	$k =~ s/^\s+//;
	$k =~ s/\s+$//;
	print "Je vais push $k\n";
	$hash_chemins{1} .= "," . $k ;
	$nb_feuilles++;

	$k = <FEUILLE> ;
	
	
	if (defined($k))
	{
		$k =~ s/^\s+//;
	$k =~ s/\s+$//;
		print "Je vais push $k\n";
		$hash_chemins{2} .= "," . $k ;
		$nb_feuilles++;
	}

	$k = <FEUILLE> ;
	
	
	if (defined($k))
	{
		$k =~ s/^\s+//;
	$k =~ s/\s+$//;
		print "Je vais push $k\n";
		$hash_chemins{3} .= "," . $k ;
		$nb_feuilles++;
	}

	$k = <FEUILLE> ;
	

	if (defined($k))
	{
		$k =~ s/^\s+//;
		$k =~ s/\s+$//;
		print "Je vais push $k\n";
		$hash_chemins{4} .= "," . $k ;
		$nb_feuilles++;
	}
	# print "k = $k\n";
}

print "toto\n";

while( my ($k,$v) = each(%hash_chemins) ) {
   print "#############Clef=$k############\n";
   my @tmp = split(',',$v);
   foreach my $t (@tmp)
   {
   	if ($t ne "")
   	{
   		print " t = $t\n";
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
# foreach my $chemin (@les_chemins)
# {
# 		$count++;
#         my $pid = fork();
#         if ($pid) {
#         # parent
#         # print "pid is $pid, parent $$\n";
#         push(@childs, $pid);
#         } elsif ($pid == 0) {
#                 # child
#                 &write($count,$chemin);
#                 exit 0;
#         } else {
#                 die "couldnt fork: $!\n";
#         }
# }


# foreach (@childs) {
#         my $tmp = waitpid($_, 0);
#          # print "done with pid $tmp\n";
# }
 
# print "End of main program\n";
########################################################
########################################################
######### END FORK #####################################
########################################################
########################################################

foreach my $chemin (@les_chemins1) {
		$count++;
        my $t = threads->new(\&write, $count,$chemin);
        push(@threads,$t);
}
foreach (@threads) {
        my $num = $_->join;
        # print "done with $num\n";
}
 
 
sub write {
        my ($num,@chemins) = (@_);

        foreach my $k (@chemins)
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
				# print "nom_fichier = $nom_fichier\n";
				# system("./../count_kmer/src/count_kmer -i $j -k $fichier_pattern -l $taille_read -o $nom_fichier -t X");
				print "./../count_kmer/src/count_kmer -i $j -k $fichier_pattern -l $taille_read -o $nom_fichier -t X >> $nom_fichier\n";
				# system("./../count_kmer/src/count_kmer -i $j -k $fichier_pattern -l $taille_read -o $nom_fichier -t X >> $nom_fichier");

			}
			close(WEKA);
		}
        sleep $num;
        # print "done with child process for $num\n";
        return $num;
}

# system("rm A954xRwm");
# print WEKA "@relation "





