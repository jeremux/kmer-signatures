#!/usr/bin/perl -w

use strict;

use Bio::DB::EUtilities;
use Getopt::Long;




###################################################################
#################### DEBUT TRAITEMENT PARAMETRE ###################
###################################################################

my $id = -1 ;
my $path = "";
my @not_id ;
my $more = "";
my $qwery = "";
my $nombre_resultat = 0;
my $mine ;
my $mail;
my $help;
my $nom_fichier;

GetOptions ('i|id=i' => \$id,
			'm|mail=s' => \$mail,
			'o|out=s'	=>\$nom_fichier,
            'path=s'   => \$path,
            'not=s'    => \@not_id,
            'more=s'	 => \$more,
            'mine'		=> \$mine,			
            'help|h' => \$help);
                   
sub usage{
    print <<HELP;

	Récupère le genbank à partir d'un taxid.
    USAGE
         perl qweryNCBI.pl     -id taxid -m youremail
			  [-path pathToSave]
			  [-not id1,id2,...]
			  [-more qwery]
			  [-mine]
			  [-help]


         -id le taxid dans la taxonomy du ncbi.
         
         -m  votre email

         -path Où sauvegarder le genbank récupéré. 

         -not Attend une liste d'id à ne pas considérer dans la requête

         -more Afin d'apporter d'autres contraintes à la requête

         -mine Considere uniquement la requête spécifiée dans more
       
         -help retourne ce message:

HELP

	exit;
}

open (LOG, '>', 'qweryNCBI.log') || die "Can't open file:$!\n";

if($help)
{
	&usage();
}

if(!$mail)
{
	print STDERR "Email obligatoire.\n";
	print STDERR "perl qweryNCBI.pl -id taxid -m email\n";
	print STDERR "perl qweryNCBI.pl -h: Pour l'aide\n";
	exit;
}

if(!$nom_fichier)
{
	print STDERR "Indiquez le nom du fichier de sortie.\n";
	print STDERR "perl qweryNCBI.pl -out name_file\n";
	print STDERR "perl qweryNCBI.pl -h: Pour l'aide\n";
	exit;
}
else
{
	$nom_fichier = $nom_fichier.'.gb';
	
	if (-e $nom_fichier)
	{
		print STDERR "Le fichier $nom_fichier existe. Specifiez un autre nom\n";
		exit;
	}	
}

if ($mine)
{ 
	if ($more ne "")
	{
		$qwery  = $more;
	}
	else
	{
		print STDERR "Precisez votre requête: perl getChild.pl -mine -more votre_requete\n";
		print STDERR "perl qweryNCBI.pl -h: Pour l'aide\n";
		exit;
	}
}
else
{
	if ($id eq -1)
	{

		print "Taxid par défaut utilisé: 2759 (Eukaryota)\n";
		$qwery = "txid2759[Organism:exp] AND (mitochondria[Title] OR mitochondrion[Title] OR mitochondrial[Title]) AND \"complete genome\"[Title]";

	
	} 
	else
	{
		$qwery = "txid".$id."[Organism:exp] AND (mitochondria[Title] OR mitochondrion[Title] OR mitochondrial[Title]) AND \"complete genome\"[Title]";
		#~ $qwery = "txid".$id."[Organism:exp] AND \"complete genome\"[Title]";
		#~ $qwery = "txid".$id."[Organism:exp]";
		
	
	}

		if ($more ne "")
		{
			$qwery = $qwery.$more;
		}

		print LOG "====================================================================================================\n";
		print LOG "=== qwery = $qwery ===\n";
		print LOG "====================================================================================================\n";
}

if (@not_id)
{
	@not_id = split /,/, join(',',@not_id);
	foreach my $x (@not_id)
	{
		#~ print "x = $x\n";
		$qwery = $qwery."(NOT txid".$x."[Organism:exp])";
	}
} 


if($path eq "")
{
		$path=".";
}
          



###################################################################
##################### FIN TRAITEMENT PARAMETRE ####################
###################################################################


###################################################################
###################### DEBUT COEUR PROGRAMME ######################
###################################################################

###################################
##### CONSTRUCTION DE L'OBJET #####
###################################

my $objetEUtil = Bio::DB::EUtilities->new(-eutil      => 'esearch',
                                       -email      => $mail,
                                       -db         => 'nucleotide',
                                       -term       => $qwery,
                                       -usehistory => 'y');
                                       
#~ my $objetEUtil = Bio::DB::EUtilities->new(-eutil      => 'esearch',
                                       #~ -email      => 'emailBidon@bidon.com',
                                       #~ -db         => 'nucleotide',
                                       #~ -term       => 'txid'.$id.'[Organism:exp] AND "complete genome"[Title]',
                                       #~ -usehistory => 'y');                                       


my $requete_correct =  $objetEUtil->get_corrected_query;

if (defined ($requete_correct))
{
	print LOG "Vouliez vous dire: \"",$requete_correct, "\"?\n";
	exit;
}

$nombre_resultat = $objetEUtil->get_count;

if (!defined($nombre_resultat))
{
	print LOG "No result for the qwery: ";
	print LOG "$qwery";
	print LOG "\n";
	exit; 
}

print LOG "Nombre de résultat attendu = $nombre_resultat\n";
print "Nombre de résultat attendu = $nombre_resultat\n";



 
# récupère l'historique
my $hist  = $objetEUtil->next_History || die 'Pas d\'historique\n';

# print "Historique retourné\n";

#########################################
##### PARAMETRE DU FORMAT DE SORTIE #####
#########################################
$objetEUtil->set_parameters(-eutil   => 'efetch', -rettype => 'fasta', -history => $hist, -verbose=> 2);

my $recommencer = 0;

##########################
##### PLUS DE LIMITE #####
##########################
my ($retmax, $retstart) = (500,0);


##########################
##### FICHIER SORTIE #####
##########################
open (my $fichier_out, '>', $path.'/'.$nom_fichier) || die "Can't open file: $!\n";

REVENIR_ICI:
	# on recupère tout
	while ($retstart < $nombre_resultat) {
	
		# remise à jour des paramètres pour contourner la limite
	    $objetEUtil->set_parameters(-retmax   => $retmax, -retstart => $retstart);	
	    # on récupère la réponse, /!\ attention à l'utilisation d'eval, fonction anonyme /!\
	    eval{
	        $objetEUtil->get_Response(-cb => sub {my ($donnes) = @_; print $fichier_out $donnes} );
	    };
	
	    # Problème lors de get_Reponse
	    if ($@) {
	        die "Erreur serveur: $@. Réessayer !\n" if $recommencer == 5;
	        print LOG "Erreur serveur, nouvel essai #$recommencer\n";
	        $recommencer++ && redo REVENIR_ICI;
	    }
	
	    $retstart += $retmax;
	}
	
close $fichier_out;

my $taille_finale = `cat $nom_fichier | grep "LOCUS" | wc -l`;
chop($taille_finale);

if ($taille_finale eq $nombre_resultat)
{
	print LOG "Toutes les séquences ont été récupérées";
	print "Toutes les séquences ont été récupérées";
}
else
{
	my $diff = $nombre_resultat - $taille_finale;
	print LOG "Warning: Il manque $diff séquences"; 
	print "Attention: Il manque $diff séquences";
}

###################################################################
####################### FIN COEUR PROGRAMME #######################
###################################################################
