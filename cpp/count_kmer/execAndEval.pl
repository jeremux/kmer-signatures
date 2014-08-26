#!/usr/bin/perl -w

use strict;

use Getopt::Long;
use Cwd 'abs_path';



my $list;
my $root;
my $start;
my $end;
my $step;
my $out;
my $nCross;
my $learn;
my $key="genomes";
my $kmerPath;
my $sample = -1;


GetOptions (
			'list=s' => \$list,
            'root=s' 	 => \$root,
            'start=i' => \$start,
            'end=i'	 => \$end,			
            'step=i' => \$step,
            'nCross=i' => \$nCross,
            'key=s' => \$key,
             'kmer=s' => \$kmerPath,
             'sample=i' => \$sample);


# exmple de ligne de commande:

# ./execAndEval.pl --root "../../create_db/Eukaryota__2759/Alveolata__33630/" --start "100" --end "300" --step "50" --nCross "10" --list list.txt --kmer "pattern.txt" --sample "20"


system("sh compilBayesJava.sh;");
system("make realclean;");
system("make;");
$root = abs_path($root);

open(LIST,'<',$list) || die "Can't open $list\n";

while (<LIST>) {
	my $l = $_;
	chomp($l);

	print "learn length = $l\n";
	if($l ne "")
	{

		
		print "./count_kmer --learn \"$l\" --sample \"$sample\" --key $key --kmer $kmerPath --root \"$root\" --start \"$start\" --step \"$step\" --end \"$end\";\n";
		system("./count_kmer --learn \"$l\" --sample \"$sample\" --key \"$key\" --kmer \"$kmerPath\" --root \"$root\" --start \"$start\" --step \"$step\" --end \"$end\";");

		if ($l eq "-1")
		{
			$l = "complete";
		}

		for (my $i = $start; $i <= $end; $i+=$step) 
		{
			system("echo  -n \"$l\t$i\t\" >> result.log;");
			print "sh execCrossVal.sh $root $nCross $i $l >> result.log;\n";
			system("sh execCrossVal.sh $root $nCross $i $l >> result.log;");
			system("echo  \"\" >> result.log;");
		}
	}
}

system("make realclean;");
