#!/usr/bin/perl

use strict;

die "perl $0 [gff file] [fasta file]\nAuthor: Yiyuan LI, liyiyuan\@genomics.org.cn
 Please cite: Xin Zhou et al, GigaScience, Ultra-deep sequencing enables high-fidelity recovery of biodiversity for bulk arthropod samples without PCR amplification.\n" unless (@ARGV == 2);

my $pos_file = shift;
my $seq_file = shift;

my %gene;

read_gff($pos_file,\%gene);

open(IN,$seq_file)||die("failed $seq_file\n");

$/=">"; <IN>; $/="\n";	
while (<IN>) {
	my $output;
	my $chr=$1 if(/^(\S+)/);
	$/=">";
	my $seq=<IN>;
	chomp $seq;
	$seq=~s/\s//g;
	$/="\n";

	my $seq_len=length($seq);
	my $chr_pp=$gene{$chr};
	foreach  my $gene (sort keys %$chr_pp) {
		my $strand=$$chr_pp{$gene}{strand};
		next if(!exists $chr_pp->{$gene}{exon});
		my @exon = @{$chr_pp->{$gene}{exon}};
		

		my $mrna;
		my ($left_leng, $right_leng);
		$left_leng = $exon[0][0]-1 if($left_leng > $exon[0][0]-1);
		$right_leng = $seq_len - $exon[-1][1] if($right_leng > $seq_len - $exon[-1][1]);

		$mrna .= substr($seq,$exon[0][0]-$left_leng-1,$left_leng) if($left_leng);
		for (my $i=0; $i<@exon; $i++) {
			$mrna .= substr($seq,$exon[$i][0]-1, $exon[$i][1] - $exon[$i][0] + 1);	
		}
		$mrna .= substr($seq,$exon[-1][1],$right_leng) if($right_leng);

		$mrna = Complement_Reverse($mrna) if($strand eq '-');
		Display_seq(\$mrna);
		my $mark = "$gene  [mRNA]  locus=$chr:$exon[0][0]:$exon[-1][1]:$strand";
		$output .= ">".$mark."\n".$mrna;
		
	}
	print $output;
}

close(IN);




####################################################
################### Sub Routines ###################
####################################################

#display a sequence in specified number on each line
#usage: disp_seq(\$string,$num_line);
#		disp_seq(\$string);
#############################################
sub Display_seq{
	my $seq_p=shift;
	my $num_line=(@_) ? shift : 50; ##set the number of charcters in each line
	my $disp;

	$$seq_p =~ s/\s//g;
	for (my $i=0; $i<length($$seq_p); $i+=$num_line) {
		$disp .= substr($$seq_p,$i,$num_line)."\n";
	}
	$$seq_p = ($disp) ?  $disp : "\n";
}
#############################################


#############################################
sub Complement_Reverse{
	my $seq=shift;
	$seq=~tr/AGCTagct/TCGAtcga/;
	$seq=reverse($seq);
	return $seq;

}
#############################################


sub read_gff{
	my $file=shift;
	my $ref=shift;
	open (IN,$file) || die ("fail open $file\n");
	while (<IN>) {
		next if(/^\#/);
		s/^\s+//;
		s/\s+$//;
		my @t = split(/\t/);
		my $tname = $t[0];
		my $qname;
		if ($t[2] eq 'mRNA' || $t[2] eq 'CDS') {
			$qname = $1 if($t[8] =~ /^ID=([^;]+);*/ || $t[8] =~ /^Parent=([^;]+);*/);
		}
		
		if ($t[2] eq 'mRNA') {
			$ref->{$tname}{$qname}{strand} = $t[6];
		}
		if ($t[2] eq 'CDS') {
			push @{$ref->{$tname}{$qname}{exon}}, [$t[3],$t[4]];
		}
	}
	close(IN);
	
}
