#!/usr/bin/perl -w

use strict;

#use Bio::DB::EUtilities;
use Getopt::Long;
use Statistics::Descriptive;

my $fichier_in = "" ;
my $titre = "" ;
my $help;

my %tableTruePositive = ();


GetOptions ('in=s'	 => \$fichier_in,
			'title=s' => \$titre,
			 'help|h' => \$help);

sub usage{
    print <<HELP;
	
USAGE
	perl make_grpah.pl -in fichier_vrais_positifs -title out_name


	-in le fichier de vrais vrais positifs avec le format
		taille_read_apprentissage	taille_read_prediction	pourcentage_vrais_positifs
	ex:
		1000	100	59.6
		1000	150	69.4
		1000	200	67.3
		1000	250	80.9
		...

	-title nom de sortie sans l'exentension 
		ex: alveolata


	USAGE 1: Genere le fichier alveolata.pdf

			perl make_grpah.pl -in result.log -title alveolata			 
HELP

	exit;
}

if($help || $fichier_in eq "" || $titre eq "")
{
	&usage();
}


open (RESULT, '<', $fichier_in) || die "Can't open file:$fichier_in\n";

my $id_learn;

while (<RESULT>) {
	
	

	my $l = $_;
	chomp($l);
	if ($l =~ /^\s*#/)
	{
		next;
	}
	else
	{
		my @columns = split(/\t/, $l);
		push(@{$tableTruePositive{$columns[0]}{$columns[1]}},$columns[2]);
	}
	
}


my $cpt = 0;
my $nb_curve = 0;
foreach my $taille_learn ( keys %tableTruePositive) 
{
	$nb_curve += 1;
	$cpt = 0;
	open (OUT, '>', "data_".$taille_learn.".txt") || die "Can't open file:$!\n";
	print OUT "#x\tbox_min\tQ1\tmedian\tq3\tmax\twidth\tlabel\n";
	foreach my $taille_predict ( sort {$a<=>$b} keys %{$tableTruePositive{$taille_learn}})
	{
		$cpt += 1;
		my @tmp_table = @{$tableTruePositive{$taille_learn}{$taille_predict}};
		{

			my $stat = Statistics::Descriptive::Full->new();
			foreach my $val (@tmp_table)
			{
				$stat->add_data($val);
			}

			my $q0 = $stat->quantile(0);
			my $q1 = $stat->quantile(1);
			my $q2 = $stat->quantile(2);
			my $q3 = $stat->quantile(3);
			my $q4 = $stat->quantile(4);
			print OUT "$cpt\t$q0\t$q1\t$q2\t$q3\t$q4\t0.3\t$taille_predict\n";
		}
	}
	close OUT;
} 


open (PLOT,'>',"boxplot.plot") || die "Can't open file:$!\n";

print PLOT "reset
set term pdf font \"Times,3\"
set output \"$titre.pdf\"
set xlabel \"read length\"   
set ylabel \"true positive\"
set bars 2.0
set style fill empty
plot ";


my $cpt2 = 0;
foreach my $taille_learn ( keys %tableTruePositive) 
{

	$cpt2 +=1 ;

	my $name_file = "data_".$taille_learn.".txt";

	if($cpt2!=$nb_curve)
	{
		if($cpt2==1)
		{
			print PLOT "'$name_file' using 1:3:2:6:5:xticlabels(8) with candlesticks title \'Quartiles $taille_learn' whiskerbars, \\
    '$name_file'         using 1:4:4:4:4 with candlesticks lt -1 notitle, \\
    '$name_file'         using 1:4 smooth bezier title 'read_size = $taille_learn', \\"
		}
		else
		{
					print PLOT "\n\t'$name_file' using 1:3:2:6:5:xticlabels(8) with candlesticks title \'Quartiles $taille_learn' whiskerbars, \\
    '$name_file'         using 1:4:4:4:4 with candlesticks lt -1 notitle, \\
    '$name_file'         using 1:4 smooth bezier title 'read_size = $taille_learn', \\"
		}

	}
	else
	{
		print PLOT "\n'$name_file' using 1:3:2:6:5:xticlabels(8) with candlesticks title \'Quartiles $taille_learn' whiskerbars, \\
    '$name_file'         using 1:4:4:4:4 with candlesticks lt -1 notitle, \\
    '$name_file'         using 1:4 smooth bezier title 'read_size = $taille_learn'"
	}
	

} 

close PLOT;

system("gnuplot boxplot.plot");

print "\n\n===============================\n";
print "==$titre.pdf generated==\n";
print "===============================\n";

