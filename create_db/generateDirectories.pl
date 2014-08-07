#!/usr/bin/perl -w

use strict;
use warnings;
use Cwd 'abs_path';

use Data::Dumper;
use Bio::DB::Taxonomy;
use Getopt::Long;

####################################################
############ Variables globales ####################
####################################################

#j'utilise deux tables de hachages pour effectuer moins d'opÃ©rations
my %tableAccession = ();
my %cptAutre = ();
my %tableOther = ();
my %tableAfficheOther = ();
my $ID_RACINE = -1;
my $path = "";
my $help;
my $time;
my $genbank = "";
my $boundary = -1;
my $flag_mkdir = 0;
my $nb_requete_entrez = 0;



GetOptions ('i|id=i' => \$ID_RACINE,
            'path=s'   => \$path,
            'time' 	 => \$time,
            'bound=i' => \$boundary,
            'gen=s'	 => \$genbank,			
            'help|h' => \$help);
            
sub usage{
    print <<HELP;
	
Génère l'arborescence taxonomique Ã  partir d'un taxon.

    USAGE
         perl generateDirectories.pl     -id taxid -bound n
					[-path pathToGenerate]
					[-time]
					[-help]
					


         -id le taxid dans la taxonomy du ncbi a 
	     partir duquel génrerer l'arborescence

         -path ou commencer cette arborescence.
	      Par defaut le dossier courant. 
         
         -bound Critere pour la creation d'un nouveau dossier
       
         -time Affiche le detail sur les temps d'execution
		 
         -help retourne ce message:

	USAGE 1: Generer l'arborescence a  partir de l'Eukaryota (taxid: 2759)
		 Avec un seuil de 20 pour la creation de nouveau dossier
			 
			 perl generateDirectories -id 2759 -bound 20
	
	USAGE 2: Genere l'arborsecencea  partir des Alveolates (taxid: 33630)
		 Avec un seuil de 10 dans /home/\$USER
			 
			 perl generateDirectories -id 33630 -bound 10 -path /home/\$USER
			 
HELP


	exit;
}

if($help)
{
	&usage();
}

if ($genbank eq "")
{
	print STDERR "Le fichier genbank est attendu:\n";
	print STDERR "perl generateDirectories -gen genomes.genbank\n";
	print STDERR "perl generateDirectories -h: Pour l'aide\n";
	exit
}

if ($ID_RACINE eq -1)
{
    print STDERR "Un taxid est attendu:\n";
	print STDERR "perl generateDirectories -id taxid";
	print STDERR "perl generateDirectories -h: Pour l'aide\n";
	exit;  
} 

if ($boundary eq -1)
{
	print STDERR "Indiquez le seuil pour la crÃ©ation d'un nouveau dossier:\n";
	print STDERR "perl generateDirectories -bound n\n";
	print STDERR "perl generateDirectories -h: Pour l'aide\n";
	exit;  
} 


# Valeur par dÃ©faut du dossier de crÃ©ation 
if ($path eq "")
{	
	$path = ".";
}


####################################################
############ BASE DONNEE LOCALE ####################
####################################################
my $dbLocale = Bio::DB::Taxonomy->new(-source   => 'flatfile',
                                  -directory=> './bdd',
                                  -nodesfile=> './bdd/nodes.dmp',
                                  -namesfile=> './bdd/names.dmp');
                                  
my $dbEntrez = Bio::DB::Taxonomy->new(-source   => 'entrez');



################################
# parametre : organism         #
# retourne : taxid de l'espece #
################################
sub get_idLocale
{
	my ($nom_espece) = @_ ;
	$nom_espece =~ s/\"//g;
	my $res = $dbLocale->get_taxonid($nom_espece);
	return $res;
}


################################
# parametre : organism         #
# retourne : taxid de l'espece #
################################
sub get_idEntrez
{
	my ($nom_espece) = @_ ;
	$nom_espece =~ s/\"//g;
	$nb_requete_entrez += 1;
	my $res = $dbEntrez->get_taxonid($nom_espece);
}


############################################
# parametre : taxid                        #
# retourne : taxid de l'ancetre            #
############################################
sub get_id_ancestorLocale
{
	my ($txid) = @_ ;
	my $taxon = $dbLocale->get_taxon(-taxonid => $txid);
	my $taxonPere = $dbLocale->ancestor($taxon);
	my $res = $taxonPere->ncbi_taxid();
	return $res;
}


sub get_id_ancestorEntrez
{
	my ($txid) = @_ ;
	my $taxon = $dbEntrez->get_taxon(-taxonid => $txid);

	# my $taxonPere2 = $dbLocale->ancestor($taxon);
	# print Dumper($taxonPere2);
	my $taxonPere = $dbEntrez->ancestor($taxon);
	my $res = $taxonPere->ncbi_taxid();
	$nb_requete_entrez += 1;
	return $res;
}




#############################################
# parametre : seoncdes                      #
# affiche: temps humainement comprÃ©hensible #
#############################################
sub FormatSeconde 
{
 
  my ($seconde) = @_;
 
  my $une_seconde = 1;
  my $une_minute  = $une_seconde * 60;
  my $une_heure   = $une_minute * 60;
  
 
  my $heure = int( $seconde / $une_heure );
  my $minute = int( ( $seconde % $une_heure ) / $une_minute );
  my $seconde2 = int( ( $seconde % $une_heure ) % $une_minute );
 

  my $retour
      = $seconde > $une_heure  ? "$heure" . " h $minute" . " min $seconde2 " . "sec"
      : $seconde > $une_minute ? $minute . " min $seconde2 " . "sec"
      : $seconde > $une_seconde ? $seconde . " secondes"
      :                           "moins d'une seconde";
  print "$retour\n";
}


##################################################################################################
##################################################################################################
################################### AFFICHAGE ####################################################
##################################################################################################
##################################################################################################

sub print_accession
{
	my ($x,$path) = @_;
	
	my @tmp_table = @{$tableAccession{$x}};

	my $cpt;

	my $taille = @tmp_table ;


	#if ($taille == $boundary)
	{
		foreach (@tmp_table)
		{

			my $k = $_;

			if(++$cpt == scalar(@tmp_table))
			{
				print GENERATE_GEN "$k";
			}
			else
			{
				print GENERATE_GEN "$k,";
			}
		} 

		print GENERATE_GEN ",$path";
		print GENERATE_GEN "\n";
	}	
}


sub print_accession_other
{
	my ($x, $racine) = @_;

	my @tmp_table = @{$tableOther{$x}};

	my $cpt;

	foreach (@tmp_table)
		{

			my $k = $_;

			if(++$cpt == scalar(@tmp_table))
			{
				print GENERATE_GEN "$k";
			}
			else
			{
				print GENERATE_GEN "$k,";
			}
		} 
}


##################################################################################################
##################################################################################################
################################### FIN AFFICHAGE ################################################
##################################################################################################
##################################################################################################

##### Ouverture fichier qui contient la liste des espÃ¨ces #####


# je grep sur ACCESSION suivi de 3 espaces pour eviter les faux accession
system "egrep -e \"ORGANISM|ACCESSION\\s\\s\\s\" $genbank > listOrganism.txt && sed -i -e \"s/\\s*ACCESSION\\s*\\|\\s*ORGANISM\\s*//g\" listOrganism.txt";
 
open (FICHIER_IN, '<', 'listOrganism.txt') || die "Can't open file:$!\n"; 

############################################
# Boucle pour traiter toutes les espÃ¨ces   #
############################################



my $time_in_traite_fichier = time;

my $tmp = 0;

my $prefix = abs_path($path);

# ICI:
while(<FICHIER_IN>) {

	# lecture accession
	my $acc = $_;

	$tmp++;

	# print "Tour $tmp\n";
	# lecture organism
	my $l = <FICHIER_IN>;

	# tab pour recupÃ©rer accession
	# il arrive que le champ accession soit de la forme
	# EF035448 DQ869278 DQ912860 DQ912861
	# On ne doit rÃ©cupÃ©rer que le premier champ EF035448
	my @tab_ligne_acc;

	# suppression des espaces et saut de ligne inutile

	chomp($l);

	# print "acc avant = $acc\n";
	chomp($acc);
	# print "acc apres = $acc\n";
	
	my $taxid_courant ;
	

	# dÃ©coupage selon les espaces du champ accession
	@tab_ligne_acc = split (/ /,$acc);

	# on recupÃ©re l'accession
	$acc = $tab_ligne_acc[0];
	# print "acc = $acc\n";
	
	# on recup le taxid_courant du champ organism
	$taxid_courant = &get_idLocale($l);
	
	## si echec en local
	if (!(defined $taxid_courant))
	{
		# redo ICI;
		print STDERR "Warning: on recupere le taxon $l sur Entrez\n";
		$taxid_courant = &get_idEntrez($l);
	}

	my @table_succes ;

	# faire un appel a get id puis tester si taxid_courant est nul
	if(defined $taxid_courant)
	{
		
		push(@table_succes,$taxid_courant);
		# $tableCompteur{$taxid_courant} += 1;
		# push(@{$tableAccession{$taxid_courant}},$acc);
		
		my $error = 0;

		# On remonte jusqu'Ã  la racine ( Eukaryota par exemple)

		
		while ( ($taxid_courant ne $ID_RACINE) && ($error eq 0))
		{	
			# il arrive que l'ancÃªtre ne soit pas retrouvÃ© en local
			eval
			{
				# on recupÃ¨re l'ancÃªtre
				my $id_ancetre = &get_id_ancestorLocale($taxid_courant);
				push(@table_succes,$id_ancetre);
				$taxid_courant = $id_ancetre;
			}
			or do
			{
				eval
				{
					print STDERR "Warning: on recupere le pere de $taxid_courant sur Entrez\n";
					my $id_ancetre = &get_id_ancestorEntrez($taxid_courant);					
					push(@table_succes,$id_ancetre);
					$taxid_courant = $id_ancetre;

				}
				# problÃ¨me en local et avec Entrez
				or do
				{
					$error = 1;
				}
			}
		}# fin while

		# aucun problÃ¨me
		if ($error eq 0)
		{
			foreach my $x (@table_succes)
			{
				# $tableCompteur{$x} += 1;
				push(@{$tableAccession{$x}},$acc);
			}
		}

		#problÃ¨me en local ou avec entrez
		else
		{

			foreach my $x (@table_succes)
			{
		
				print STDERR "Trace erreur: Get ancestor of the taxon with taxid $x: KO\n";
			}
		}
	}

	else
	{
		print STDERR "WARNING: problÃ¨me pour rÃ©cupÃ©rer l'ID de $l avec accession = $acc\n";
	}
	
	
}


if ($time)
{
	$time_in_traite_fichier = time - $time_in_traite_fichier;
	print "Temps traitement comptage : ";
	&FormatSeconde($time_in_traite_fichier);
}	

	
##### Fichier script_mkdir.sh #####
open (SCRIPT_MKDIR, '>', 'script_mkdir.sh') || die "Can't open file:$!\n";
open (SCRIPT_MKDIR_DATA, '>', 'script_mkdir_data.sh') || die "Can't open file:$!\n";

##### Fichier generateGenbank.sh #####
open (GENERATE_GEN,  '>', 'listGenbank2.txt') || die "Can't open file:$!\n"; 


# En tete des scripts
print SCRIPT_MKDIR "#!/bin/sh\n";





# on recupere l'objet taxon
my $objet_taxon = $dbLocale->get_taxon(-taxonid => $ID_RACINE);



my $flag_first = 0;
my $flag_tmp = 1;
my $flag_clean = 0;

my $path_racine = "";
my $first = 0;

sub traite_taxon
{
	# on recupere le taxon passÃ© en paramÃ¨tre
	my ($taxon_courant) = @_ ;

	my $flag_other = 0;
	my $flag_un = 0;

	my $flag_enfant = 0;

	# on recupÃ¨re les enfants
	my @lesEnfants = $dbLocale->each_Descendent($taxon_courant);

	my $nom_taxon = $taxon_courant->scientific_name;
	$nom_taxon =~ s/ /_/g;
	my $taxid = $taxon_courant->ncbi_taxid();


	
	if ($first == 0)
	{
		$path_racine =  $path."/".$nom_taxon."__$taxid";
		$path = $prefix . "/" . $path."/".$nom_taxon."__$taxid";
		$first = 1;
	}
	else
	{
		$path =  $path."/".$nom_taxon."__$taxid";
	}
	
	$path =~ s/\(//g;
	$path =~ s/\)//g;
	$path =~ s/\s+/_/g;


	# faire en dehors de la fonction 
	
	my $racine = $path;

	print SCRIPT_MKDIR "mkdir -p ";
	print SCRIPT_MKDIR "$path\n";

	#######################################
	######### Classer les données #########
	#######################################
	print SCRIPT_MKDIR_DATA "mkdir -p ";
	print SCRIPT_MKDIR_DATA "$path/data/fasta/nucleotides\n";
	print SCRIPT_MKDIR_DATA "mkdir -p ";
	print SCRIPT_MKDIR_DATA "$path/data/fasta/aminoAcids\n";
	print SCRIPT_MKDIR_DATA "mkdir -p ";
	print SCRIPT_MKDIR_DATA "$path/data/genbank/\n";
	print SCRIPT_MKDIR_DATA "mkdir -p ";
	print SCRIPT_MKDIR_DATA "$path/frequencies/\n";

	# print SCRIPT_MKDIR "/$nom_taxon\n";

	# print GENERATE_GEN "perl extractGenbank.pl -gen $genbank -list ";
	&print_accession($taxid,$path);
	#print GENERATE_GEN ",$path";
	#print GENERATE_GEN "\n";
	
	foreach my $enfant (@lesEnfants)
	{
		my $taxIdEnfant = $enfant->ncbi_taxid();

		my $nb_genome = $#{$tableAccession{$taxIdEnfant}}+1;

		if (defined $nb_genome && $nb_genome >= $boundary)
		{
			# $path = $path."/".$nom_taxon;
			&traite_taxon($enfant);
			$flag_enfant = 1;
			$path = $racine;
		}
		else
		{
			if (defined $nb_genome )
			{
				if ($nb_genome > 0)
				{
					if ($flag_other eq 0)
					{
						$flag_other = 1;
					}
					
					push(@{$tableOther{$taxid}},@{$tableAccession{$taxIdEnfant}}); 
				}
			}
			else
			{
				if($flag_un eq 0)
				{
					# &print_accession($taxid,$path);
					print "$path\n";
					$flag_un = 1;
				}
			}
		}
	}

	# si on a crÃ©Ã© un dossier ($flag_enfant eq 1)
	# et si on possÃ¨de un other ($flag_other eq 1)
	if (($flag_other eq 1) && ($flag_enfant eq 1))
	{
		print SCRIPT_MKDIR "mkdir -p ";
		print SCRIPT_MKDIR "$racine";
		print SCRIPT_MKDIR "/others\n";

		#######################################
		######### Classer les données #########
		#######################################
		print SCRIPT_MKDIR_DATA "mkdir -p ";
		print SCRIPT_MKDIR_DATA "$racine";
		print SCRIPT_MKDIR_DATA "/others/data/fasta/nucleotides\n";
		print SCRIPT_MKDIR_DATA "mkdir -p ";
		print SCRIPT_MKDIR_DATA "$racine";
		print SCRIPT_MKDIR_DATA "/others/data/fasta/aminoAcids\n";
		print SCRIPT_MKDIR_DATA "mkdir -p ";
		print SCRIPT_MKDIR_DATA "$racine";
		print SCRIPT_MKDIR_DATA "/others/data/genbank/\n";
		print SCRIPT_MKDIR_DATA "mkdir -p ";
		print SCRIPT_MKDIR_DATA "$racine";
		print SCRIPT_MKDIR_DATA "/others/frequencies/\n";
	

					
		# print GENERATE_GEN "perl extractGenbank.pl -gen $genbank -list ";
		&print_accession_other($taxid,$racine);
		print GENERATE_GEN ",$racine/others";
		print GENERATE_GEN "\n";
	}
}
	

my $cpt_autre = 0;
my $profondeur = 0;
my $flag_ecrit = 0;

sub tree
{
	my ($taxon_courant) = @_ ;

	$profondeur++;

	# my $flag_other = 0;

	my $flag_enfant = 0;

	# on recupÃ¨re les enfants
	my @lesEnfants = $dbLocale->each_Descendent($taxon_courant);

	my $nom_taxon = $taxon_courant->scientific_name;

	$nom_taxon =~ s/\(//g;
	$nom_taxon =~ s/\)//g;

	my $taxid = $taxon_courant->ncbi_taxid();

	my $nb_genome_pere = $#{$tableAccession{$taxid}}+1;
	print KRONA "<node name=\"$nom_taxon\">\n
					<genomes><val>$nb_genome_pere</val></genomes>\n";

	
	if ($flag_ecrit eq 0)
	{
		print NEWICK "(";
	}
	else
	{
		print NEWICK ",(";
		$flag_ecrit = 0;
	}

	
	foreach my $enfant (@lesEnfants)
	{
		$cpt_autre = 0;
		

		my $taxIdEnfant = $enfant->ncbi_taxid();

		my $nom_enfant = $enfant->scientific_name;

		my $nb_genome = $#{$tableAccession{$taxIdEnfant}}+1;

		if (defined $nb_genome && $nb_genome >= $boundary)
		{
			&tree($enfant);
			$profondeur--;
			$flag_enfant = 1;
			
		}
		else
		{
			if (defined $nb_genome && $nb_genome > 0)
			{
				$cpt_autre+=$nb_genome;
				# print "hello\n";
				# $flag_other = 1;
			}
		}

		
	}

	$cpt_autre = $#{$tableOther{$taxid}}+1;
	# si on a crÃ©Ã© un dossier ($flag_enfant eq 1)
	# et si on possÃ¨de un other ($flag_other eq 1)
	if (($cpt_autre gt 0) && ($flag_enfant eq 1))
	{
		print NEWICK ",others[";
		print NEWICK "$cpt_autre";
		print NEWICK "]:";
		$profondeur++;
		print NEWICK "$profondeur";
		$profondeur--;
		# print NEWICK ")";
		print KRONA "<node name=\"others\">\n
					<genomes><val>$cpt_autre</val></genomes>\n
					</node>";
		
	}

	print KRONA "</node>";
	print NEWICK ")";
	print NEWICK "$nom_taxon";
	print NEWICK "[";
	print NEWICK "$nb_genome_pere";
	print NEWICK "]:";
	print NEWICK "$profondeur";
	$flag_ecrit = 1;

}


# print NEWICK "(";

my $nom_racine = $objet_taxon->scientific_name;
$nom_racine =~ s/ /_/g;
my $nom_fichier = '../generate_data/generateGenbank_'.$nom_racine.'.sh';
open (SCRIPT_GET, '>',$nom_fichier) || die "Can't open file:$!\n";

open (NEWICK, '>', 'tree_'.$nom_racine.'.newick') || die "Can't open file:$!\n"; 


print SCRIPT_GET "#!/bin/sh\n";

open (KRONA, '>', 'krona_'.$nom_racine.'.xml') || die "Can't open file krona.xml!\n";

print KRONA "<krona>\n
<attributes magnitude= \"genomes\">\n
<attribute display=\"genomes\">genomes</attribute>\n
</attributes>\n
	<color attribute=\"genomes\" valueStart=\"0\" valueEnd=\"55\" hueStart=\"120\" hueEnd=\"240\"></color>\n
	<datasets>\n
		<dataset>Brand X</dataset>
	</datasets>";


if ($time)
{
	my $time_in_traite_taxon = time;
	&traite_taxon($objet_taxon);
	&tree($objet_taxon);
	$time_in_traite_taxon = time - $time_in_traite_taxon;

	print "Temps pour generer les scripts: ";
	&FormatSeconde($time_in_traite_taxon);
}
else
{
	&traite_taxon($objet_taxon);
	&tree($objet_taxon);
}

my $abs_path = abs_path($path_racine);
my @les_chemins = split (/\//,$abs_path);


my $cpt2;

my $taille = @les_chemins ;

$abs_path = "";

foreach (@les_chemins)
{

	my $k = $_;

	if(++$cpt2 != scalar(@les_chemins))
	{
		$abs_path .= "$k\/"
	}
}

$abs_path =~ s/ /_/g;
print "abs_path = $abs_path\n";

my $abs_path_racine = $abs_path . $path_racine;
$abs_path_racine =~ s/ /_/g;

my $abs_genbank = abs_path($genbank);
$abs_genbank =~ s/ /_/g;
print SCRIPT_GET "perl extractGenbank.pl -list listGenbank.txt -gen $abs_genbank -conf conf --root $abs_path;\n";
print SCRIPT_GET "bash fillAll_v2.sh ";
print SCRIPT_GET "2> 46Ukt6xMyK6Li6 .txt ;";
print SCRIPT_GET "rm -rf 46Ukt6xMyK6Li6.txt;";

# print "path_racine = $abs_path\n";

print NEWICK ";";
system "sed -i -e s/\\(\\)//g tree.newick";
print "Nombre de requete sur Entrez : $nb_requete_entrez\n";

print KRONA "</krona>";
close KRONA;

close SCRIPT_MKDIR;
close GENERATE_GEN;
close SCRIPT_MKDIR_DATA;
system "chmod +x script_mkdir.sh";
system "chmod +x script_mkdir_data.sh";
system "chmod +x $nom_fichier";
system "rm -rf $abs_path_racine && ./script_mkdir.sh && ./get_leaf.sh $abs_path_racine listGenbank2.txt $nom_racine && ./script_mkdir_data.sh && rm -f listOrganism.txt";

# system "./script_mkdir_data.sh"
#system "rm -rf listGenbank2.txt";

# &affiche();
