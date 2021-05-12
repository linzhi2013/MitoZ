#! usr/bin/perl -w
=head1 Copyright

raw_reads_filter_v0.5.pl

Copyright (c) 2017-2018 Guanliang Meng <mengguanliang@foxmail.com>.

This file is part of MitoZ.

MitoZ is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MitoZ is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MitoZ.  If not, see <http://www.gnu.org/licenses/>.

=cut

use strict;
use Getopt::Long;

=head1 Description

	Usage: perl raw_reads_filter.pl <parameter>

	-1 <STR>       input read1 fastq file
	-3 <STR>       output read1 fastq file
	-k <BEG,END>   keep only bases between BEG and END for each read [full length read]
	-n <INT>       the maximun N allowed in single end reads [10]
	-q <QUAL,CONT> the maixmum percentage of low quality (ASCII code's integer)
	               bases allowed in single end reads [50,40]
	-z             output gzip result
	-h             print out this information

	Filtering steps:
	    1. trim reads (-k);
	    2. filter Ns reads (-n);
	    3. filter low quality reads (-q).

	Version: 0.5 for SE
	Update: use md5 value of read1-read2 as hash key to decrease the Mem cost
	Author: mengguanliang@genomics.cn

=cut

my ($fq1, $fq3, $keep_range, $help, $nma, $qua_cont, $zip);

## there is a trap: if we assign $nma eq 0, then `$nma || =10;`  will be excuted, the $nma get 10 now!!
$nma = 10;
$qua_cont = "50,40";

GetOptions(
	"1=s"=>\$fq1,
	"3=s"=>\$fq3,
	"k:s"=>\$keep_range,
	"n:i"=>\$nma,
	"q:s"=>\$qua_cont,
	"z!"=>\$zip,
	"h|help"=>\$help
);
die `pod2text $0` if ($help|| !$fq1 || !$fq3);

my ($low_qua, $cont) = split /,/, $qua_cont;

print("input files:\n\t$fq1\n");
print("output files:\n\t$fq3\n");
print("filter parameters:\n");
print("\tk, keep only bases between BEG and END for each read: ", ($keep_range)?($keep_range):('full length read'), "\n");
print("\tn, the maximun N allowed in single end reads: $nma\n");
print("\tq, the maixmum percentage of low quality (ASCII code's integer) bases allowed in single end reads: $qua_cont\n");

my $rd_len = 0;

my $is_N = 0;
my $is_lowqua = 0;
my $clean_rd = 0;


###################

# open file handles
#open(FIN, "<gzip(autopop)", $fq1) || die $!;
#open(RIN, "<gzip(autopop)", $fq2) || die $!;
if($fq1=~/\.gz$/){
	open(FIN, "gzip -dc $fq1 |") || die $!;
}else{
	open(FIN, "<", $fq1) || die $!;
}

if($zip){
	if($fq3!~/.gz$/){
		$fq3 .= ".gz";
	}
	open(F_O, "| gzip > $fq3") || die $!;
}else{
	open(F_O, ">", $fq3) || die $!;
}

my ($f01, $f02, $f03, $f04);
my ($k_start, $k_end, $k_len);

if($keep_range){
		($k_start, $k_end) = split /,/, $keep_range;
		$k_len = $k_end - $k_start + 1;
		$k_start = $k_start - 1;
}

while($f01=<FIN>){
	chomp($f01);
	chomp($f02=<FIN>);
	chomp($f03=<FIN>);
	chomp($f04=<FIN>);


	if($keep_range){
		$f02 = substr($f02, $k_start, $k_len);
		$f04 = substr($f04, $k_start, $k_len);
	}


	if($f02=~tr/Nn/Nn/ > $nma){
		$is_N++;
		next;
	}

	if(&filter_qua($f04, $low_qua, $cont)){
		$is_lowqua++;
		next;
	}

	$clean_rd++;
	print F_O "$f01\n$f02\n$f03\n$f04\n"; #print "$f01\n" if ($clean_rd == 1);
}

close FIN;
close F_O;


###################################

###################################

sub filter_qua{
	my ($qua, $low_qua, $cont) = @_;
	my $count = 0;
	my $q;
	my $seq_len = length($qua);
	for (my $i=0; $i<=$seq_len-1; $i++){
		$q = substr($qua, $i, 1);
		$count++ if (ord($q) < $low_qua);
	}
	if ($count*100 / $seq_len > $cont){
		return 1;
	}else{
		return 0;
	}
}

