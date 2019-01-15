#! /usr/local/bin/perl -w
use strict;
# declaration variables
our ($seqfile, $com, $nosolar, $solar_bin, $optional);
my ($FileIn, $FileOut);
die "Usage: perl $0 [input -m 8 format]\n Author: Yiyuan LI, liyiyuan\@genomics.org.cn  
 Please cite: Xin Zhou et al, GigaScience, Ultra-deep sequencing enables high-fidelity recovery of biodiversity for bulk arthropod samples without PCR amplification.\n" unless (@ARGV == 1);
# initialization
my $tmp = `dirname $0`;
chomp $tmp;
$solar_bin = "$tmp/solar";

my $opt_n = 100000;
my $opt_d = -1;
my $opt_c = 1;
my $opt_C = 1;

$com = "| $solar_bin";
$com .= " -n $opt_n";
$com .= " -d $opt_d";
$com .= " -c";
$com .= " -C";

open (FILEIN, $ARGV[0]) || die "fail $ARGV[0]";
open(FILEOUT, $com);
$FileOut = \*FILEOUT;

while (<FILEIN>) {
        my @t = split;
        print $FileOut $t[0],"\t0\t",$t[1],"\t0\t",$t[6],"\t",$t[7],"\t",$t[8],"\t",$t[9],"\t";
	print $FileOut int($t[11] + 0.5),"\t",$t[10],"\n";
}

close FILEIN;
