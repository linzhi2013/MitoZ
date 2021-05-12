#!/usr/bin/perl -w
=head1 Copyright

extract_no-mapping-reads_v2.pl

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

To extract mapped or no-mapped reads according to mapped SAM file (which contains only
mapped reads). gzip supported.

updates:

When extracting no-mapped reads (custom output):

  1. This script starts to work only when the max_length of sequences
     in sam file >= Min_len cutoff

  2. Only extract unmapped reads for the max_len_seq

  3. Also extract mapped reads on the 5'/3' ends of the max_len_seq


Usage: perl extract_no-mapping-reads.pl <parameter>

  -h, --help  show this help message and exit
  -sam <str>  mapping SAM file which contains only the mapping reads
  -q1 <str>   in fastq 1 file
  -q2 <str>   in fastq 2 file
  -q3 <str>   out fastq 1 file
  -q4 <str>   out fastq 2 file
  -map        output all mapped reads
  -unmap      output all unmapped reads
  -z          gzip out? default: False

for custom output:

  -rl <int>   read length (default: 150)
  -m  <int>   minimum max_len of sequences (default: 10000)
  -a  <int>   the length of 5' ends to be extracted mapped reads (default: 170)
  -b  <int>   the length of 3' ends to be extracted mapped reads (default: 170)
=cut

my ($sam, $q1, $q2, $q3, $q4, $map, $unmap, $read_len, $zip, $min_len, $region_5, $region_3, $help);

$read_len = 150;
$min_len = 10000;
$region_5 = 170;
$region_3 = 170;

GetOptions(
	"sam=s" => \$sam,
	"q1=s"  => \$q1,
	"q2=s"  => \$q2,
	"q3=s"  => \$q3,
	"q4=s"  => \$q4,
	"-map!" => \$map,
	"-unmap!" => \$unmap,
	"rl:i"   => \$read_len,
	"m:i"   => \$min_len,
	"a:i"   => \$region_5,
	"b:i"   => \$region_3,
	"z!"    => \$zip,
	"h|help" => \$help,
);

die `pod2text $0` if ($help || !$sam || !$q1 || !$q2 || !$q3 || !$q4);

my (%seqids, $seqid, $scaf_name, $max_len, $scaf_kept, $seq_start);
$max_len = 0;

$min_len = 0;

open(IN, "<", $sam) || die $!;
while (<IN>) {
	if(/^\@/){
		if (/\@SQ\s+SN:(.+?)\s+LN:(\d+)/){
			if ($2 > $max_len){
				$scaf_kept = $1;
				$max_len = $2;
			}
		}
		next;
	}

	if ($map || $unmap){
		$seqids{$seqid} = 1;
		next;
	}else{
		if ($max_len < $min_len){
			exit(0);
		}
	}

	($seqid, $scaf_name, $seq_start) = (split /\s+/)[0, 2, 3];

	if ($scaf_name eq $scaf_kept){
		if (0 < $seq_start and ($seq_start + $read_len) < $region_5){
			next;
		}elsif (($max_len - $region_3) < $seq_start){
	 		next;
	 	}else{
			$seqids{$seqid} = 1;
		}

	}
}
close IN;

if($q1=~/\.gz$/){
	open (FIN, "gzip -dc $q1 |") || die $!;
	open (RIN, "gzip -dc $q2 |") || die $!;
}else{
	open (FIN, "<", $q1) || die $!;
	open (RIN, "<", $q2) || die $!;
}

if($zip){
	if($q3 !~ /.gz$/){
		$q3 .= ".gz";
		$q4 .= ".gz";
	}
	open(F_O, "| gzip > $q3") || die $!;
	open(R_O, "| gzip > $q4") || die $!;
}else{
	open(F_O, ">", $q3) || die $!;
	open(R_O, ">", $q4) || die $!;
}

my ($f01, $f02, $f03, $f04, $r01, $r02, $r03, $r04);
while ($f01=<FIN>){
	$f02=<FIN>;
	$f03=<FIN>;
	$f04=<FIN>;

	$r01=<RIN>;
	$r02=<RIN>;
	$r03=<RIN>;
	$r04=<RIN>;

	$seqid = (split /\s+/, $f01)[0];
	$seqid =~ s/^@//;

	if ($map and exists($seqids{$seqid})){
		print F_O "$f01$f02$f03$f04";
		print R_O "$r01$r02$r03$r04";
		next;

	}elsif($unmap and !exists($seqids{$seqid})){
		print F_O "$f01$f02$f03$f04";
		print R_O "$r01$r02$r03$r04";
		next;
	}

	unless (exists($seqids{$seqid})){
		print F_O "$f01$f02$f03$f04";
		print R_O "$r01$r02$r03$r04";
	}

}

close F_O;
close R_O;
