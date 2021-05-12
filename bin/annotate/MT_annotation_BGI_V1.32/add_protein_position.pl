#!usr/bin/perl -w
use strict;

die "Usage: perl $0 [.solar.filter.nr] [.cds]\n" unless (@ARGV == 2);
open (NR, $ARGV[0]) or die "$ARGV[0] $!\n";
my %hash;
while(<NR>){
	chomp;
	my @line = split /\s+/;
	$hash{$line[0]} = $line[0]."_Ref".$line[2].":".$line[3]."aa";	
}


open (FA, $ARGV[1]) or die "$ARGV[1] $!\n";
open OUT, ">$ARGV[1].position.cds" or die "$ARGV[1].position.cds $!\n";
while(<FA>){
       	chomp;
	if(/^>/){
		my @name_line = split;
		$name_line[0] =~ s/\>//;
		print OUT ">$hash{$name_line[0]}\t$name_line[1]\t$name_line[2]\n";
	}else{
		print OUT "$_\n";
	}
}

close NR;
close FA;
close OUT;
print "DONE!";
