#!/usr/bin/perl -w

use strict;

use Bio::DB::EUtilities;
use Getopt::Long;

my $help;
my $time;
my $conf;
my $genbank = "";
my $generate = "";
my $list_acc ="";
my $boundary = -1;
my $flag_mkdir = 0;
my $path =".";
my $racine = "";

my %hash_genbank = ();
my %hash_name = ();
my %hash_sequence = ();
my %hash_gene = ();
my %hash_alias = ();
my %hash_protein = ();


GetOptions ('gen=s'	 => \$genbank,
            'list=s' => \$generate,
            'conf=s' => \$conf,			
            'help|h' => \$help,
            'root=s' => \$racine);

# if($racine eq "")
# {
# 	print STDERR "La racine de la bdd locale est attendu:\n";
# 	print STDERR "perl generateDirectories -root path_to_racine\n";
# 	print STDERR "perl generateDirectories -h: Pour l'aide\n";
# 	exit
# }
# else
# {
# 	my $copy_racine = $racine;
# 	$copy_racine = substr $copy_racine , -1 ; #dernier caractère
# 	print "copy_racine = $copy_racine\n";
# 	if ($copy_racine ne "/")
# 	{
# 		$racine = $racine."/";
# 	}
# }

# fichier d'entré

open (GENBANK, '<', $genbank) || die "Can't open file:$!\n";
open (LIST, '<', $generate) || die "Can't open file:$!\n";
open (ALIAS, '<', $conf ) || die "Can't open file:$!\n";
open (LOG, '>', 'extractGenbank.log') || die "Can't open file:$!\n";
open (LOG2, '>', 'extractGenbank.log2') || die "Can't open file:$!\n";

# fichier généré

while (<ALIAS>)
{
	my $l = $_;
	chomp($l);
	my @ligne = split (/:/,$l);
	my $nom_gene = $ligne[0];
	$nom_gene =~ s/^\s+|\s+$//g;

	$nom_gene = lc $nom_gene;
	$hash_alias{$nom_gene} = $nom_gene;

	# $ligne[1] =~ s/^\s+|\s+$//g;

	if (defined($ligne[1]))
	{
		my @alias = split(/,/,$ligne[1]);
		foreach my $x (@alias)
		{
			$x =~ s/^\s+|\s+$//g;
			$x = lc $x;
			$hash_alias{$x} = $nom_gene;
			
		}
	}
}


my $fin = "//";


my $case = "";
my $seq = "";
my $grep_accession = "\\^ACCESSION\\s+";
my $grep_definition = "DEFINITION\\s+";
my $grep_organism = "ORGANISM";
my $grep_origin = "ORIGIN";
my $grep_cds = "\\s+CDS\\s+";
my $grep_protein = "translation=";
my $grep_complement = "complement";
my $grep_guillemet ="\"";
my $flag_origin = 0;
my $flag_take_protein = 0;
my $flag_take_gene = 0;
my $gene;
my $protein = "";

# my $grep_accession = "LOCUS\\s\\s\\s\\s\\s\\s\\s";

my $acc;
my $organism;
my @tab_lesProtein;

# decoupage en bloc

my $flag_definition = 0;
my $flag_cds = 0;
my $flag_complement = 0;
my $position ;
my $nom_standard;

my $cpt = 0;
my $cpt_ligne = 0; 
my $taille_gen = 0;

my $taille_finale = `cat $genbank | grep "LOCUS" | wc -l`;

while(<GENBANK>)
{
	my $l = $_;

	$cpt_ligne++;


	my $copy_l = $l;
	chomp($l);
	# $l =~ s/^\s+|\s+$//g;

	# if ($cpt_ligne > 8818000 && $cpt_ligne < 8824895)
	# {
	# 	print "$l\n";
	# }
	# if ($cpt_ligne > 8824895 )
	# {
	# 	exit;
	# }


	if($l =~ m/^\/\/$/)
	{
		# print DEBUG "===========\n";
		# print DEBUG "$fin\n";
		$cpt++;
		$flag_origin = 0;
		chomp($case);
		chomp($acc);
		$acc =~ s/^\s+|\s+$//g;
		# $acc =~ s/\s+/_/g;
		
		print LOG2 "$acc\n";
	
		
		$hash_genbank{$acc} = $case;
		$taille_gen = $taille_gen + 1;
		$hash_name{$acc} = $organism . "__" . $acc;
		$hash_sequence{$acc} = $seq;


		$case = "";
		$seq = "";

		print "\033[2J";
		print "\033[0;0H";
		my $rap = ($cpt / $taille_finale)*100;
		my $rap_string = "" . $rap;
		$rap_string = substr($rap_string,0,4);
		print "Traitement 1 / 3\n";
		print "Avancement: $rap_string %\n";

	}
	else
	{
		# print DEBUG "$l\n";
		$case = $case . $l . "\n";

		# $l =~ s/^\s+|\s+$//g;

		if ($flag_origin eq "1")
		{

			# trim
			
			$l =~ s/\d//g;
			$l =~ s/\s*//g;


			$seq = $seq . $l . "\n";
		}
		else
		{
			if ($l =~ m/(\s\sORGANISM\s*)(.*)/ )
			{
				
				$organism = $2;

				# trim
				$organism =~ s/^\s+|\s+$//g;

				$organism =~ s/\s+/_/g;

				$organism =~ s/\(//g;
				$organism =~ s/\)//g;
				
				
			}

			if ($l =~ m/$grep_origin/)
			{
				$flag_origin = 1;
			}

			##########################
			##########################

			if ($l =~ m/$grep_definition/)
			{
				# $case = $case . $l . "\n";
				$flag_definition = 1;
			}

			if ($l =~ m/(ACCESSION\s+)(.*)/ && $flag_definition eq "1")
			{
				# print "l = $l\n";
				# $l =~ s/$sed_accession//g;
				# print "acc = $l\n";
				# $acc = $l;
				########################
				########################

				$acc = $2;
				
				

				my @tab_ligne_acc;

				# découpage selon les espaces du champ accession
				@tab_ligne_acc = split (/ /,$acc);

				# on recupére l'accession
				$acc = $tab_ligne_acc[0];

				# print "acc = $acc\n";

				$flag_definition = 0; 
				$acc =~ s/^\s+|\s+$//g;
				# chop($acc);
				

			}



			if ($l =~ m/$grep_cds/)
			{
				$flag_cds = 1;

				my $ligne_cds = "";
				my $ligne = "" ;

				$l =~ s/^\s+|\s+$//g;
				$ligne = $l;
				$position = "";

				# ##########################
				# print "ligne avant = $ligne\n";
				while(!($ligne =~ m/(\/)(.*)/))
				{

					
					chomp($ligne);
					$position = $position . " " . $ligne;

					my $second_lect = <GENBANK>;
					$cpt_ligne++;
				
					$case = $case . $second_lect . "\n";

					chomp($second_lect);

					# 	if ($cpt_ligne > 8818000 && $cpt_ligne < 8824895)
					# {
					# 	print "$second_lect\n";
					# }
					# if ($cpt_ligne > 8824895 )
					# {
					# 	exit;
					# }	
					
					
					$second_lect =~ s/^\s+|\s+$//g;
					$ligne = $second_lect;
					# print "ligne boucle = $l\n";

					
				}

				############
				### TODO ###

				$l = $ligne;
				# print "ligne boucle = $l\n";
				# print "position = $position\n";

				##########################
				# while(!($l =~ m/(\/gene=")(.*)(")/) || !($l =~ m/(\/note=")(.*)(")/))
				while(!($l =~ m/(\/gene=")(.*)(")/))
				# while(!($l =~ m/(\/gene=")(.*)/) || !(($l =~ m/(\/note=")/) || !($l =~ m/(\/product=")(.*)(")/))
				{

					if($l =~ m/(\/note=")/)
					{
						last;
					}
					else
					{
						if ($l =~ m/(\/product=")(.*)(")/)
						{
							last;
						}
					}
					# print "ligne = $l\n";
					$l =~ s/^\s+|\s+$//g;
					$ligne = $ligne . $l;
					my $second_lect = <GENBANK>;
					$cpt_ligne++;	
					$l = $second_lect;

					$case = $case . $second_lect . "\n";

					chomp($second_lect);

					# if ($cpt_ligne > 8818000 && $cpt_ligne < 8824895)
					# {
					# 	print DEBUG "$cpt_ligne: $second_lect\n";
					# }
					# if ($cpt_ligne > 8824895 )
					# {
					# 	exit;
					# }	
					
				}

				
				# print "arret = $l\n";

				if ($l =~ m/(\/gene=")(.*)(")/)
				{
					$gene = $2;
					# print "gene = $gene\n";
				}
				elsif ($l =~ m/(\/note=")(.*)/)
				{	
					# $gene =~ s/\"//g;
					$gene = $2;
				}
				elsif ($l =~ m/(\/product=")(.*)(")/)
				{
					$gene = $2;
				}
				else
				{
					print STDERR "Error gene non encountered\n";
				}

				$gene = lc $gene;
				$gene =~ s/\"//g;
				# print "gene = $gene\n";
				
				if (defined($hash_alias{$gene}))
				{
					$nom_standard = $hash_alias{$gene};
					$flag_take_gene = 1;
					
				}
				else
				{
					$flag_take_gene = 0;
					print LOG "WARNING: gene '$gene' avec accession $acc non reconnu/considéré\n";
				}


				



			}


			if ($flag_take_protein eq 1)
			{
				$l =~ s/^\s+|\s+$//g;
				if ($l =~ m/$grep_guillemet/)
				{
					# print "$l\n";
					$l =~ s/"//g;

					$protein = $protein . $l;
					$flag_take_protein = 0;

					$position =~ s/\s*CDS\s+//g;
					$position =~ s/^\s+|\s+$//g;

					# print "position = $position\n";
					if ($position =~ m/^complement\(/)
					{
						$position =~ s/^complement\(|\)$//g;
						$position .= "_" . "1";

					}
					elsif ($position =~ m/^join\(/)
					{
						$position =~ s/^join\(|\)$//g;
						$position .= "_" . "2";
						
					}
					else
					{	
						$position =~ s/^\(|\)$//g;
						$position .= "_" . "0";
					}

					# print "POSITION = $position\n";

					$protein = $protein . "_" . $nom_standard . "_" . $position;


					# print "debut ajout\n";
					# print "acc = $acc\n";
					push(@{$hash_protein{$acc}},$protein);
					# print "fin ajout\n";



				}
				else
				{
					$protein = $protein . $l;
				}
			}

			if (($l =~ m/$grep_protein/) && ($flag_take_gene eq 1))
			{
				$l =~ s/^\s+|\s+$//g;
				$l =~ s/"/-/g ;
				
				
				$flag_take_protein = 1;
				$protein = "";
				$protein = `echo -n $l | cut -d'-' -f2`;


				chomp($protein);
				# print "protein pour $acc = $protein \n\n\n";
			}
		
		}
		
	}
}

sub print_accession
{
	my ($x) = @_;
	
	# print "x = $x";
	my @tmp_table = @{$hash_protein{$x}};

	my $cpt;

	foreach (@tmp_table)
	{
		my $k = $_;

		print "Protein = $k\n\n";
	}
}

sub print_cut
{
	my ($dna) = @_;

	chomp($dna);


    $dna =~ s/^\s+|\s+$//g;
	
	my $res = "";
	my @tab = split(//,$dna);
	my $last = scalar(@tab) - 1;
	my $cpt = -1;
	my $k;

	for ($k = 0; $k <= $last; $k++)
	{
		if ($tab[$k] ne ' ' && $tab[$k] ne '')
		{
			if ($cpt > 80 )
			{
				$res = $res . "\n" . $tab[$k];
				$cpt=0;
			}
			else
			{
				$cpt++;
				$res = $res . $tab[$k];
				
			}
		}
	}

	return $res; 

}


my $flag_affiche = 0 ;

my @tab_extract = <LIST>;



# Pour chaque morceau x du fichier genbank

my $path_fasta = $path; 
my $path_protein;
my $path_nt;

sub reverse_complement_IUPAC 
{
	my ($dna) = @_;

	chomp($dna);
	$dna =~ s/^\s+|\s+$//g;
	# reverse the DNA sequence
	my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
	$revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
	return $revcomp;
}

sub reverse_complement 
{
    my ($dna) = @_;


    my $c;
    chomp($dna);
    $dna =~ s/^\s+|\s+$//g;
    # print "on va prendre le complement de: \n";
    # print "$dna";
    my $res = "" ;
    my $k;

    my @tab = split(//,$dna);
    my $last = scalar(@tab) - 1;
    
    # http://www.bioinformatics.org/sms/iupac.html
    for ($k = $last; $k >= 0; $k--)
    {
    	$c = $tab[$k];
    	$c = lc $c;

    	if($c eq "a") 
    	{
  			$res = $res . "t";
		} elsif($c eq "c") {
  			$res = $res . "g";
		} elsif($c eq "g") {
  			$res = $res . "c";
		} elsif($c eq "t" || $c eq "u") {
  			$res = $res . "a";
  		} elsif($c eq "n") {
  			$res = $res . "n";
  		} elsif($c eq "y") {
  			$res = $res . "r";
		} elsif($c eq "r") {
  			$res = $res . "y";
		} elsif($c eq "k") {
  			$res = $res . "m";
		} elsif($c eq "m") {
  			$res = $res . "k";
		} else {
		    print LOG "Error in sequence : $dna while trying to reverse-complement: $c doesn\'t recognized\n\n" ;
		}
	}
    	
    return $res;
}

sub traite_normal
{
	my ($s,$pos) = @_;

	# print "=========\n";
	# 		print "normal = $pos\n";
	# 		print "=========\n";
		#on enlève les sauts de lignes et les blancs
	$s =~ s/\n//g;
	$s =~ s/ //g;

	$position =~ s/^\s+|\s+$//g;
	$pos =~ s/>|<//g;
	my @les_pos = split(/\.\./,$pos);

	my $debut = $les_pos[0];
	my $fin = $les_pos[1];
	my $long = int($fin) - int($debut) + 1;

	$debut -=1;
	my $res = substr($s,$debut,$long);
	return $res;
}

sub traite_complement
{
	my ($s,$pos) = @_;

	my @les_pos = split(/,/,$pos);

	my $res = "";
		#on enlève les sauts de lignes et les blancs
	$s =~ s/\n//g;
	$s =~ s/ //g;
	foreach my $x (@les_pos)
	{

		$x =~ s/^\s+|\s+$//g;

		# print "=========\n";
		# print "complement = $x\n";
		# print "=========\n";

		if ($x =~ m/^complement\(/)
		{

			# print "=========\n";
			# print "complement = $x\n";
			# print "=========\n";

			$res =~ s/^\s+|\s+$//g;
			$x =~ s/^complement\(|\)$//g;
			$res = $res . &traite_complement($s,$x);


			# print "position = $position\n";
		}
		elsif ($x =~ m/^join/)
		{
			# print "la ligne = $ligne\n";

			# print "=========\n";
			# print "join = $x\n";
			# print "=========\n";
			$x =~ s/^\s+|\s+$//g;
			$x =~ s/^join\(|\)$//g;
			$res = $res . &traite_join($s,$x);

			# print "position = $position\n";
		}
		else
		{	
			# print "=========\n";
			# print "normal complement = $x\n";
			# print "=========\n";
			$x =~ s/>|<//g;
			$x =~ s/^\(|\)$//g;
			my @indice = split(/\.\./,$x);
			my $debut = $indice[0];
			my $fin = $indice[1];
			my $long = int($fin) - int($debut) + 1;
			$debut = $debut-1;

	

			$res = $res . substr($s,$debut,$long);	

		}
	}

	$res = reverse_complement_IUPAC($res);

	return $res;
}


sub traite_join
{
	my ($s,$pos) = @_;

	# print "pos = $pos\n";
	# print "s = $s\n";
	my @les_pos = split(/,/,$pos);
	#on enlève les sauts de lignes et les blancs
	$s =~ s/\n//g;
	$s =~ s/ //g;
	# print "Traitement de $pos\n";

	my $res = "";
	foreach my $x (@les_pos)
	{
		# print "x = $x\n";
		$x =~ s/^\s+|\s+$//g;

		if ($x =~ m/^complement\(/)
		{
			$x =~ s/^\s+|\s+$//g;
			$x =~ s/^complement\(|\)$//g;
			$res = $res . &traite_complement($s,$x);

			# print "position = $position\n";
		}
		elsif ($x =~ m/^join/)
		{
			# print "la ligne = $ligne\n";
			
			$x =~ s/^\s+|\s+$//g;
			$x =~ s/^join\(|\)$//g;

			
			$res = $res . &traite_join($s,$x);

			# print "position = $position\n";
		}
		else
		{
			# print "Traitement de x = $x: \n";
			$x =~ s/>|<//g;
			my @indice = split(/\.\./,$x);
			
			
			my $debut = $indice[0];
			my $fin = $indice[1];

			if (!defined($fin))
			{
				$fin = $debut;
			}
			my $long = int($fin) - int($debut) + 1;

			$debut -= 1;
			# print "res avant = $res\n";
			$res = $res . substr($s,$debut,$long);	

			
			# print "res après = $res\n";
		}
	}

	return $res;
}

open (FEUILLE,  '>', 'feuille.txt') || die "Can't open file feuille.txt:$!\n";
print FEUILLE "$racine\n";

my $cpt2;
my $cpt1;
my $id;
my $id2;
my $cpt3 = 0;



foreach my $un_genbank (@tab_extract)
{

	$cpt3 = $cpt3 + 1;
	print "\033[2J";
	print "\033[0;0H";
	my $rap = ($cpt3 / $taille_gen)*100;
	my $rap_string = "" . $rap;
	$rap_string = substr($rap_string,0,4);
	print "Traitement 2 / 3 (Génération des donnees aux feuilles)\n";
	print "Avancement: $rap_string %\n";

	chomp($un_genbank);
	my $i = 0;
	my @tab_acc =  split (/,/,$un_genbank);

	# print "un_genbank = $un_genbank\n\n\n\n";
	
	my $dernier_indice ;
	my $new_seq ;

	$dernier_indice = scalar(@tab_acc) - 1;

	# $path = $racine.$tab_acc[$dernier_indice];
	$path = $tab_acc[$dernier_indice];

	print FEUILLE "$path\n";
	# print "racine = $racine \n";
	$path_fasta = $path;
	$path_protein = $path;
	$path_nt = $path;
		
	### DEBUT V2
	my @les_chemins = split(/\//,$path);
	$cpt2 = 0;
	$cpt1 = 0;
	foreach my $bout (@les_chemins)
	{
		if(++$cpt1 == scalar(@les_chemins)-1)
		{
			if ($bout =~ m/(.*)(__)(.*)/)
			{
				$id2 = $3
			}
		}
		if(++$cpt2 == scalar(@les_chemins))
		{
			# print "bout = $bout\n";
			if($bout eq "others")
			{
				$id = $id2 . "-others";
			}
			else
			{
				if ($bout =~ m/(.*)(__)(.*)/)
				{
					$id = $3
				}
			}
		}
	}
	### FIN V2
	# print "id = $id\n";

	#######################################
	######### CLASSER LES DONNEES #########
	#######################################
	$path = $path . "/data/genbank/" ;
	system "mkdir -p $path";
	$path = $path . "genomes". "_$id" . ".genbank";
	$path_fasta = $path_fasta . "/data/fasta/nucleotides/genomes/";
	system "mkdir -p $path_fasta";
 	$path_fasta = $path_fasta . "genomes" . "_$id" . ".fasta"; 
	$path_protein = $path_protein . "/data/fasta/aminoAcids/" ;
	$path_nt = $path_nt . "/data/fasta/nucleotides/";
	#######################################
	#######################################
	#######################################
	# $path = $path . "/" . "genomes". "_$id" . ".genbank";
	# $path_fasta = $path_fasta . "/" . "genomes" . "_$id" . ".fasta"; 
	# $path_protein = $path_protein . "/" ;
	# $path_nt = $path_nt . "/";


	my $copy_aa = $path_protein;
	my $copy_nt = $path_nt;

	# print "On va ouvrir $path\n";
	open (GENERATE_GEN,  '>', $path) || die "Can't open file $path:$!\n"; 
	open (FASTA,  '>', $path_fasta) || die "Can't open file $path_fasta:$!\n"; 

	# Pour chaque accession tmp de la liste d'accession
	for ($i = 0; $i < $dernier_indice; $i++) 
	{	
		my $access = $tab_acc[$i];
		chomp($access);
		$access =~ s/^\s+|\s+$//g;
		# $access =~ s/\s+/_/g;
		# print "Accession = $access\n\n";
		
		if (!defined($hash_genbank{$access}))
		{
			print LOG "Error to get access: $access\n";
			next;
		}
		my $tmp = $hash_genbank{$access};
		my $name = $hash_name{$access};
		my $la_seq = $hash_sequence{$access};
		my $copie_seq = $la_seq;
		
		###############################################
		#### Si un jour on veut couper la séquence ####
		###############################################


		# my $seq_taille = length($la_seq);
		# if ($seq_taille > 140000)
		# {
		# 	$la_seq = substr $la_seq, 0, 140000;
		# }
		
		# my $cut = &print_cut($la_seq);
		print GENERATE_GEN "$tmp";
		print GENERATE_GEN "\n$fin\n\n";

		# print "Accession = $access\n\n";

		print FASTA ">";
		print FASTA "$name\n";
		# print FASTA "$cut\n\n";
		print FASTA "$la_seq\n\n";

		$la_seq =~ s/\n//g;

		if (defined(@{$hash_protein{$access}}))
		{
			my @tmp_table = @{$hash_protein{$access}};

			foreach my $aa (@tmp_table)
			{
				# print "Accession = $access\n";
				# print "aa = $aa\n\n";
				my @aa_nom = split(/_/,$aa);
				my $nom_fichier = $aa_nom[1];
				my $la_proteine = $aa_nom[0];
			
				my $take_complement = $aa_nom[3];

			
				
				####################
				#  AA              #
				####################
				#######################################
				######### CLASSER LES DONNEES #########
				#######################################
				system "mkdir -p $path_protein/$nom_fichier/";
				$path_protein = $path_protein. "/". $nom_fichier. "/" . "aa_" . $nom_fichier . "_$id" . ".fasta";
				#######################################
				#######################################
				#######################################
				# $path_protein = $path_protein . "aa_" . $nom_fichier . "_$id" . ".fasta";

				$la_proteine = &print_cut($la_proteine);
				open (PROT,  '>>', $path_protein) || die "Can't open file $path_protein:$!\n"; 
				print PROT ">";
				print PROT "aa_";
				print PROT "$name\n";
				print PROT "$la_proteine\n\n";



				# $debut -= 1;

				####################
				#  NT              #
				####################
				#system "mkdir -p $path_nt/$nom_fichier";
				#######################################
				######### CLASSER LES DONNEES #########
				#######################################
				system "mkdir -p $path_nt/$nom_fichier/";
				$path_nt = $path_nt . "/". $nom_fichier. "/". "nt_" . $nom_fichier . "_$id" . ".fasta";
				#######################################
				#######################################
				#######################################

				# $path_nt = $path_nt . "nt_" . $nom_fichier . "_$id" . ".fasta";
				open (NT,  '>>', $path_nt) || die "Can't open file $path_protein:$!\n";
				print NT ">";
				print NT "nt_";
				print NT "$name\n";

				
				if ($take_complement eq "0")
				{
					$new_seq = &traite_normal($copie_seq,$aa_nom[2]);
				}
				elsif ($take_complement eq "1")
				{
					$new_seq = &traite_complement($copie_seq,$aa_nom[2]);
				}
				else
				{
					$new_seq = &traite_join($copie_seq,$aa_nom[2]);
				}

				$new_seq = &print_cut($new_seq);
				print NT "$new_seq\n\n";
				# print PROT "$sous_sequence\n\n";
		

				$path_protein = $copy_aa;
				$path_nt = $copy_nt;
			}
		}

		
		
	}

	close GENERATE_GEN;
}


