#!/usr/bin/perl 
use strict;

die "Usage:$0 <genewise> <pep_len>\n Author: Yiyuan LI, liyiyuan\@genomics.org.cn
 Please cite: Xin Zhou et al, GigaScience, Ultra-deep sequencing enables high-fidelity recovery of biodiversity for bulk arthropod samples without PCR amplification.\n" if @ARGV<1;

my $genewise=shift;
my $pep=shift;

my %Len;
open IN,$pep or die "$!";
while(<IN>){
	if (/(\S+)\s+(\S+)/){
		$Len{$1}=$2;
	}
}
close IN;

open IN,$genewise or die "$!";
$/="//\n";
while(<IN>){
	next if !/^Bits/;
	chomp;
	my $line1=$_;
	my $line2=<IN>;
	chomp $line2;
	my $line3=<IN>;	
	chomp $line3;
	my @l1=split(/\n/,$line1);
	my $len;
	my $pid;
	my $id;
	if ($l1[1]=~/\S+\s+(\S+)(-D\d+)\s+(\d+)\s+(\d+)/){
		$id=$1.$2;
		$pid=$1;
		$len=abs($4-$3)+1;
	}elsif($l1[1]=~/\S+\s+(\S+)\s+(\d+)\s+(\d+)/){
		$id=$1;
		$pid=$1;
		$len=abs($3-$2)+1;
	}else{
		die ;
	}
	my $cover=sprintf("%.2f",$len/($Len{$pid}*100));
	my @l3=split(/\n/,$line3);
	my $shift=-1;
	my @cds;
	my $chr;
	foreach (@l3){
		my @c=split(/\t/);
		my $start;
		if ($c[0]=~/^(\S+)_(\d+)_(\d+)/){
			$chr=$1;
			$start=$2;
		}else{
			die;
		}
		@c[3,4] =@c[4,3] if $c[3] >$c[4];
		if ($c[2] eq 'match'){
			$shift++;	
		}elsif ($c[2] eq 'cds'){
			$c[3]+=$start-1;
			$c[4]+=$start-1;
			push @cds,[$c[3],$c[4],$c[6],$c[7]];
		}
	}
	@cds=sort {$a->[0]<=>$b->[0]} @cds;
	my @mRNA=($chr,'GeneWise','mRNA',$cds[0][0],$cds[-1][1],$cover,$cds[0][2],'.',"ID=$id;Shift=$shift;");
	print join("\t",@mRNA)."\n";
	foreach (@cds){
		print join("\t",$chr,'GeneWise','CDS',$$_[0],$$_[1],'.',$$_[2],$$_[3],"Parent=$id;")."\n";
	}
}
close IN;
