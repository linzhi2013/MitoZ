#!/usr/bin/perl

use strict;

die "perl $0 [cds fas]\n For the MT_annotation_BGI.pl\n Author: Yiyuan LI, liyiyuan\@genomics.org.cn\n Please cite: Xin Zhou et al, GigaScience, Ultra-deep sequencing enables high-fidelity recovery of biodiversity for bulk arthropod samples without PCR amplification.\n" unless (@ARGV == 1);
my %CODE = (
				'GCA' => 'A', 'GCC' => 'A', 'GCG' => 'A', 'GCT' => 'A',                               # Alanine
				'TGC' => 'C', 'TGT' => 'C',                                                           # Cysteine
				'GAC' => 'D', 'GAT' => 'D',                                                           # Aspartic Acid
				'GAA' => 'E', 'GAG' => 'E',                                                           # Glutamic Acid
				'TTC' => 'F', 'TTT' => 'F',                                                           # Phenylalanine
				'GGA' => 'G', 'GGC' => 'G', 'GGG' => 'G', 'GGT' => 'G',                               # Glycine
				'CAC' => 'H', 'CAT' => 'H',                                                           # Histidine
				'ATC' => 'I', 'ATT' => 'I',                                             # Isoleucine
				'AAA' => 'K', 'AAG' => 'K',                                                           # Lysine
				'CTA' => 'L', 'CTC' => 'L', 'CTG' => 'L', 'CTT' => 'L', 'TTA' => 'L', 'TTG' => 'L',   # Leucine
				'ATG' => 'M', 'ATA' => 'M',                                                           # Methionine
				'AAC' => 'N', 'AAT' => 'N',                                                           # Asparagine
				'CCA' => 'P', 'CCC' => 'P', 'CCG' => 'P', 'CCT' => 'P',                               # Proline
				'CAA' => 'Q', 'CAG' => 'Q',                                                           # Glutamine
				'CGA' => 'R', 'CGC' => 'R', 'CGG' => 'R', 'CGT' => 'R', 'AGG' => 'R',   # Arginine
				'TCA' => 'S', 'TCC' => 'S', 'TCG' => 'S', 'TCT' => 'S', 'AGC' => 'S', 'AGT' => 'S', 'AGA' => 'S', 'AGG' => 'S', # Serine
				'ACA' => 'T', 'ACC' => 'T', 'ACG' => 'T', 'ACT' => 'T',                               # Threonine
				'GTA' => 'V', 'GTC' => 'V', 'GTG' => 'V', 'GTT' => 'V',                               # Valine
				'TGG' => 'W', 'TGA' => 'W',                                                                        # Tryptophan
				'TAC' => 'Y', 'TAT' => 'Y',                                                           # Tyrosine
				'TAA' => 'U', 'TAG' => 'U',                                                           # Stop
	);


$/=">"; <>; $/="\n";
while (<>) {
	my $head = $_;
	chomp $head;
	my $key = $1 if($head =~ /^(\S+)/);
	my $phase = ($head =~ /\s+phase[:\s]+([012])\s+/i) ? $1 : 0 ;
	$/=">";
	my $seq = <>;
	chomp $seq;
	$/="\n";
	
	my $prot = cds2aa($seq,$phase);
	Display_seq(\$prot);
	print ">$head [translate_table: mitochondrion]\n".$prot;

}
close IN;


####################################################
################### Sub Routines ###################
####################################################


#display a sequence in specified number on each line
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


## translate CDS to pep
####################################################
sub cds2aa {
	my $seq = shift;
	my $phase = shift || 0;	

	$seq =~ s/\s//g;
	$seq = uc($seq);
	
	my $len = length($seq);
	
	my $prot;
	for (my $i=$phase; $i<$len; $i+=3) {
		my $codon = substr($seq,$i,3);
		last if(length($codon) < 3);
		$prot .= (exists $CODE{$codon}) ? $CODE{$codon} : 'X';
	}
	$prot =~ s/U$//;
	return $prot;

}

