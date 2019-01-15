#!/usr/bin/perl

use strict;
die "Usage: perl $0 [bait file] [database file]\n Author: Yiyuan LI, liyiyuan\@genomics.org.cn
 Please cite: Xin Zhou et al, GigaScience, Ultra-deep sequencing enables high-fidelity recovery of biodiversity for bulk arthropod samples without PCR amplification.\n" unless (@ARGV == 2);

my $BaitSeperator = '\s+'; ## default seperator of bait file
my $FishSeperator = '\s+'; ## default seperator of fish file
my ($Baitformat,$FishFormat,$BaitColumn,$FishColumn,$Except,$PatternMode,$GeneMode);
my ($Verbose,$Help);

$BaitColumn = 1;
$FishColumn = 1;

my $bait_file=$ARGV[0];
my $fish_file=$ARGV[1];

my %Bait;


read_table($bait_file,\%Bait,$BaitColumn,$BaitSeperator);

##print Dumper \%Bait;
output_table($fish_file,\%Bait,$FishColumn,$FishSeperator);


####################################################
################### Sub Routines ###################
####################################################


sub read_table{
	my ($file,$bait_hp,$bait_colum,$bait_seperator) = @_;
	open(IN,$file)||die("fail to open $file\n");
	while (<IN>) {
		chomp;
		s/^\s+//g;
		my @temp=split(/$bait_seperator/,$_);
		my $id = $temp[$bait_colum-1];

		if (!$GeneMode) {
			$$bait_hp{$id}=1;
		}else{
			$id = $1 if($id =~ /^(\w+)-\w+$/);
			$$bait_hp{$id} = 1;
			
		}
		
	}
	close(IN);
}


sub output_table{
	my ($file,$fish_hp,$fish_colum,$fish_seperator) = @_;
	open(IN,$file)||die("fail to open $file\n");
	while (<IN>) {
		chomp;
		s/^\s+//g;
		my @temp=split(/$fish_seperator/,$_);
		my $id = $temp[$fish_colum-1];
		$id = $1 if($GeneMode && $id =~ /^(\w+)-\w+$/);

		if (!$Except) {
			print $_."\n" if (!$PatternMode && exists $$fish_hp{$id});
			print $_."\n" if ($PatternMode && word_pattern_hash($id,$fish_hp));

		}else{
			print $_."\n" if (!$PatternMode && !exists $$fish_hp{$id});
			print $_."\n" if ($PatternMode && !word_pattern_hash($id,$fish_hp));
		}
		
	}
	close(IN);
}


sub word_pattern_hash{
	my $word = shift;
	my $hash_p = shift;
	my $value = 0;
	foreach my $hash_key (keys %$hash_p) {
		if ($word =~ /$hash_key/){
			$value = 1;
			last;
		}
	}
	return $value;
}
