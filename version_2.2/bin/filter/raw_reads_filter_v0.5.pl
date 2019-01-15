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
	-2 <STR>       input read2 fastq file
	-a <STR>       input 1.adapter.list.gz file
	-b <STR>       input 2.adapter.list.gz file
	-3 <STR>       output read1 fastq file
	-4 <STR>       output read2 fastq file
	-d             filter duplications (caused by PCRs After adapter ligation),
	               thus, directions of read1 & read2 of duplications should be completely identical
	-m <INT>       cut-off adapter mis-match bases [3]
	-l <INT>       cut-off adapter align length [15]
	-k <BEG,END>   keep only bases between BEG and END for each read [full length read]
	-n <INT>       the maximun N allowed in single end reads [10]
	-q <QUAL,CONT> the maixmum percentage of low quality (ASCII code's integer)
	               bases allowed in single end reads [50,40]
	-z             output gzip result
	-h             print out this information

	Filtering steps:
	    1. filter duplications (-d);
	    2. trim reads (-k);
	    3. filter adapater contamination (-a, -b, -m, -l. Not to filter if "-a" and "-b" are not specified);
	    4. filter Ns reads (-n);
	    5. filter low quality reads (-q).
	thus, if we choose trimming reads, then some reads in `1.adapter.list.gz` and
	`2.adapter.list.gz` should not be discarded if the adapter alignment
	positions are out of range (BEG, END), BUT this program discard them too
	for simplification.

	Version: 0.5
	Update: use md5 value of read1-read2 as hash key to decrease the Mem cost
	Author: mengguanliang@genomics.cn

=cut

my ($fq1, $fq2, $adapter1, $adapter2, $fq3, $fq4, $dup, $keep_range, $help, $nma, $qua_cont, $mis, $aln, $zip);

## there is a trap: if we assign $nma eq 0, then `$nma || =10;`  will be excuted, the $nma get 10 now!!
$nma = 10;
$qua_cont = "50,40";
$mis = 3;
$aln = 15;

GetOptions(
	"1=s"=>\$fq1,
	"2=s"=>\$fq2,
	"a=s"=>\$adapter1,
	"b=s"=>\$adapter2,
	"3=s"=>\$fq3,
	"4=s"=>\$fq4,
	"d!"=>\$dup,
	"m:i"=>\$mis,
	"l:i"=>\$aln,
	"k:s"=>\$keep_range,
	"n:i"=>\$nma,
	"q:s"=>\$qua_cont,
	"z!"=>\$zip,
	"h|help"=>\$help
);
die `pod2text $0` if ($help|| !$fq1 || !$fq2 || !$fq3 || !$fq4);

my ($low_qua, $cont) = split /,/, $qua_cont;

print("input files:\n\t$fq1\n\t$fq2\n");
print("\t$adapter1\n\t$adapter2\n") if ($adapter1 and $adapter2);
print("output files:\n\t$fq3\n\t$fq4\n");
print("filter parameters:\n");
print("\tk, keep only bases between BEG and END for each read: ", ($keep_range)?($keep_range):('full length read'), "\n");
print("\tm, cut-off adapter mis-match bases: $mis\n") if ($adapter1 and $adapter2);
print("\tl, cut-off adapter align length: $aln\n") if ($adapter1 and $adapter2);
print("\tn, the maximun N allowed in single end reads: $nma\n");
print("\tq, the maixmum percentage of low quality (ASCII code's integer) bases allowed in single end reads: $qua_cont\n");
if(defined($dup)){
	print("\td, Filtering duplications: True\n\n");
}else{
	print("\td, Filtering duplications: False\n\n");
}

my (%no_need_reads);

my $rd_len = 0;

my $tot_rd = 0;
my $is_N = 0;
my $is_adapter = 0;
my $is_lowqua = 0;
my $clean_rd = 0;

## filter duplication
my ($fq3_dup, $fq4_dup);
my $dup_rd = 'NA';

if(defined($dup)){
	$fq3_dup = $fq3.".dup.tmp";
	$fq4_dup = $fq4.".dup.tmp";
	($tot_rd, $dup_rd) = &filter_dup($fq1, $fq2, $fq3_dup, $fq4_dup, 0);
	$fq1 = $fq3_dup;
	$fq2 = $fq4_dup;
}

###################

# filter adapter
if ($adapter1 and $adapter2){
	&filter_adapter($adapter1, \%no_need_reads);
	&filter_adapter($adapter2, \%no_need_reads);
}

# open file handles
#open(FIN, "<gzip(autopop)", $fq1) || die $!;
#open(RIN, "<gzip(autopop)", $fq2) || die $!;
if($fq1=~/\.gz$/){
	open(FIN, "gzip -dc $fq1 |") || die $!;
	open(RIN, "gzip -dc $fq2 |") || die $!;
}else{
	open(FIN, "<", $fq1) || die $!;
	open(RIN, "<", $fq2) || die $!;
}

if($zip){
	if($fq3!~/.gz$/){
		$fq3 .= ".gz";
		$fq4 .= ".gz";
	}
	open(F_O, "| gzip > $fq3") || die $!;
	open(R_O, "| gzip > $fq4") || die $!;
}else{
	open(F_O, ">", $fq3) || die $!;
	open(R_O, ">", $fq4) || die $!;
}

my ($f01, $f02, $f03, $f04, $r01, $r02, $r03, $r04, $f_out, $r_out);
my ($k_start, $k_end, $k_len);
my ($f01_id, $r01_id);

if($keep_range){
		($k_start, $k_end) = split /,/, $keep_range;
		$k_len = $k_end - $k_start + 1;
		$k_start = $k_start - 1;
}

while($f01=<FIN>){
	unless (defined($dup)){
		$tot_rd++;
	}
	chomp($f01);
	chomp($f02=<FIN>);
	chomp($f03=<FIN>);
	chomp($f04=<FIN>);

	chomp($r01=<RIN>);
	chomp($r02=<RIN>);
	chomp($r03=<RIN>);
	chomp($r04=<RIN>);

	if($keep_range){
		$f02 = substr($f02, $k_start, $k_len);
		$f04 = substr($f04, $k_start, $k_len);
		$r02 = substr($r02, $k_start, $k_len);
		$r04 = substr($r04, $k_start, $k_len);
	}

	#$f01_id = (split /\t/, $f01)[0];
	#$r01_id = (split /\t/, $r01)[0];
	if((exists $no_need_reads{$f01}) || (exists $no_need_reads{$r01})){
		$is_adapter++;
		next;
	}
	if($f02=~tr/Nn/Nn/ > $nma || $r02=~tr/Nn/Nn/ > $nma){
		$is_N++;
		next;
	}
	if(&filter_qua($f04, $low_qua, $cont) || &filter_qua($r04, $low_qua, $cont)){
		$is_lowqua++;
		next;
	}

	$clean_rd++;
	print F_O "$f01\n$f02\n$f03\n$f04\n"; #print "$f01\n" if ($clean_rd == 1);
	print R_O "$r01\n$r02\n$r03\n$r04\n"; #print "$r01\n" if ($clean_rd == 1);;
}

close FIN;
close RIN;
close F_O;
close R_O;

$rd_len = length($f02);

print("#raw_read_bases\tduplication\tis_adapter\ttoo_many_Ns\tis_low_quality\tclean_read_bases\tclean_read_length\n");
if($dup_rd eq 'NA'){
	print($tot_rd*$rd_len*2,"\t", $dup_rd,"\t",$is_adapter*$rd_len*2,"\t", $is_N*$rd_len*2, "\t", $is_lowqua*$rd_len*2, "\t", $clean_rd*$rd_len*2, "\t", $rd_len, "\n");
}else{
	print($tot_rd*$rd_len*2,"\t", $dup_rd*$rd_len*2,"\t",$is_adapter*$rd_len*2,"\t", $is_N*$rd_len*2, "\t", $is_lowqua*$rd_len*2, "\t", $clean_rd*$rd_len*2, "\t", $rd_len, "\n");
}

if(defined($dup)){
	unlink($fq3_dup);
	unlink($fq4_dup);
}

###################################

###################################

sub filter_adapter{
	my ($adapter_list_file, $ref_no_need_reads) = @_;
	my ($read_id, $align_len, $mismatch);

	if ($adapter_list_file=~/\.gz$/){
		open IN, "gzip -dc $adapter_list_file |" or die $!;
	}else{
		open IN, "<", $adapter_list_file or die $!;
	}
	while(<IN>){
		next if(/^#/);
		my @line = split /\t/;
		($read_id, $align_len, $mismatch) = @line[0, -2, -1];
		$read_id = '@'.$read_id; # reads id in 1.adapter_list_file do not have '@' character.
		$$ref_no_need_reads{$read_id} = 1 if ($align_len >= $aln && $mismatch <= $mis);
	}
	close IN;
	return $ref_no_need_reads;
}

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

sub filter_dup{
	#use PerlIO::gzip;
	use Digest::MD5 qw(md5);
	# in theory, we can use other better algorithms, e.g. Digest::SHA, but I think MD5 is just fine here.

	my ($fq1, $fq2, $fq3, $fq4, $zip) = @_;
	my %hash;
	my $tot_rd = 0;
	my $out_rd = 0;
	my $dup_rd = 0;
	my $dup_rate = 0.0;
	my ($f01, $f02, $f03, $f04, $r01, $r02, $r03, $r04);
	my $md5 = "";
	#my ($fr, $fr_r, $fr_c, $fr_cr, $r02_r);

	#open(FIN, "<gzip(autopop)", $fq1) || die $!; ## fail to open large gzip file, always stop on 185th read.
	#open(RIN, "<gzip(autopop)", $fq2) || die $!;
	if($fq1=~/\.gz$/){
		open(FIN, "gzip -dc $fq1 |") || die $!;
		open(RIN, "gzip -dc $fq2 |") || die $!;
	}else{
		open(FIN, "<", $fq1) || die $!;
		open(RIN, "<", $fq2) || die $!;
	}
	if($zip){
		if($fq3!~/.gz$/){
			$fq3 .= ".gz";
			$fq4 .= ".gz";
		}
		open(F_O, "|gzip>$fq3") || die $!;
		open(R_O, "|gzip>$fq4") || die $!;
	}else{
		open(F_O, ">", $fq3) || die $!;
		open(R_O, ">", $fq4) || die $!;
	}
	while($f01=<FIN>){
		$tot_rd++;# print("tot_rd: $tot_rd\n");
		chomp($f01);
		chomp($f02=<FIN>);
		chomp($f03=<FIN>);
		chomp($f04=<FIN>);

		chomp($r01=<RIN>);
		chomp($r02=<RIN>);
		chomp($r03=<RIN>);
		chomp($r04=<RIN>);

		$md5 = md5("$f02$r02");
		unless (exists $hash{$md5}){
			print F_O "$f01\n$f02\n$f03\n$f04\n";
			print R_O "$r01\n$r02\n$r03\n$r04\n";
			$hash{$md5} = 1;
		}

		#$fr = $f02.$r02;
		# filter duplications (caused by PCRs After adapter ligation),
        # thus, directions of read1 & read2 of duplications should be completely identical
		#$r02_r = reverse($r02);
		#$fr = $f02.$r02_r;
		#$fr_r = reverse($fr);

		#$fr_c = $fr; $fr_c =~tr/ATGCNatgcn/TACGNtacgn/;
		#$fr_cr = reverse($fr_c);

		#unless (exists $hash{$fr} && exists $hash{$fr_r} && exists $hash{$fr_c} && exists $hash{$fr_cr}){
		#	print F_O "$f01\n$f02\n$f03\n$f04\n";
		#	print R_O "$r01\n$r02\n$r03\n$r04\n";
		#	$hash{$fr} = 1;
		#}
	}

	close FIN;
	close RIN;
	close F_O;
	close R_O;

	$out_rd = keys %hash;
	$dup_rd = $tot_rd - $out_rd;

	#$dup_rate = $dup_rd*100 / $tot_rd;
	#$dup_rate = sprintf("%.2f%%", $dup_rate);
	#print("#raw_read_pairs\tduplication_read_pairs\tduplication_rate\toutput_read_pairs\n");
	#print("$tot_rd\t$dup_rd\t$dup_rate\t$out_rd\n\n");

	return ($tot_rd, $dup_rd);
}

