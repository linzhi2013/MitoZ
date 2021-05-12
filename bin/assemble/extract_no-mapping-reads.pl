#!/usr/bin/perl -w
=head1 Copyright

extract_no-mapping-reads.pl

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

To extract no mapping reads according to mapped SAM file (which contains only
mapping reads). gzip supported.


Usage: perl extract_no-mapping-reads.pl <parameter>

  -h, --help  show this help message and exit
  -sam <str>  mapping SAM file which contains only the mapping reads.
  -q1 <str>   in fastq 1 file.
  -q2 <str>   in fastq 2 file.
  -q3 <str>   out fastq 1 file.
  -q4 <str>   out fastq 2 file.
  -z          gzip out? default: False
=cut

my ($sam, $q1, $q2, $q3, $q4, $zip, $help);

GetOptions(
	"sam=s" => \$sam,
	"q1=s"  => \$q1,
	"q2=s"  => \$q2,
	"q3=s"  => \$q3,
	"q4=s"  => \$q4,
	"z!"    => \$zip,
	"h|help"=> \$help,
);

die `pod2text $0` if ($help || !$sam || !$q1 || !$q2 || !$q3 || !$q4);

my (%seqids, $seqid);

open(IN, "<", $sam) || die $!;
while (<IN>) {
	next if(/^@/);
	$seqid = (split /\s+/)[0];
	$seqids{$seqid} = 1
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

	unless (exists($seqids{$seqid})){
		print F_O "$f01$f02$f03$f04";
		print R_O "$r01$r02$r03$r04";
	}

}

close F_O;
close R_O;
