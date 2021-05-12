#!/usr/bin/perl -w
use strict;

die "usage: perl $0 [ contig file ]\nAuthor: Yiyuan LI, liyiyuan\@genomics.org.cn  
 Please cite: Xin Zhou et al, GigaScience, Ultra-deep sequencing enables high-fidelity recovery of biodiversity for bulk arthropod samples without PCR amplification.\n" unless (@ARGV==1);
open IN,$ARGV[0] or die "$!\n";
my($head,$seq);
my @line;
my $line;
$/=">";<IN>;$/="\n";
my (@length, @seq_length);

while ($line=<IN>) {
	@line=split(/\s+/,$line);
        my $head=$line[0];
        chomp $head;        
        $/=">";
        my $seq=<IN>;
        chomp $seq;
	$seq=~s/\s+//g;
	$seq=~s/>//g;  
      	my $len=length($seq);

     	print "$head\t$len\n";
        $/="\n";
}
close IN;
