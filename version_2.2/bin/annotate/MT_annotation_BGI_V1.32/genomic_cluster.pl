#!/usr/bin/perl

use strict;

die "perl $0 [input file filter.table]\n Author: Yiyuan LI, liyiyuan\@genomics.org.cn
 Please cite: Xin Zhou et al, GigaScience, Ultra-deep sequencing enables high-fidelity recovery of biodiversity for bulk arthropod samples without PCR amplification.\n" unless (@ARGV==1);

my $Overlap_percent = 0.5;
my $InFile=shift;
my %genes;

#table:
#At1g01010.1     Chr12   +       1138400 1146949  cds_len    score
#$t[0]	        $t[5]    $t[4]   $t[7]   $t[8]    $len       $t[10]

Read_Table($InFile,\%genes);

foreach my $chr(sort keys %genes) {
	my @sort_genes= sort {$a->[3]<=>$b->[3]} @{$genes{$chr}};
	my @cluster;
	foreach my $gene (@sort_genes) {
		if (@cluster==0) {
			push @cluster,[@$gene];
		}else{
			my $overlap;
			$overlap=find_overlap(\@cluster,$gene,$Overlap_percent,'percent');
			if ($overlap==1) {
				push @cluster,[@$gene];
			}else{
				my @noredundance;
				find_nonredundance(\@cluster, \@noredundance);
				@cluster=();
				push @cluster,[@$gene];

				foreach my $gene (@noredundance){
					print join("\t",@$gene)."\n";
				}
			}
		}
	}
	my @noredundance;
	find_nonredundance(\@cluster, \@noredundance);
	@cluster=();
	foreach my $gene (@noredundance){
		print join("\t", @$gene)."\n";
	}
}




#At1g01010.1     Chr12   +       1138400 1146949
sub Read_Table {
	my $InFile=shift;
	my $genes_p=shift;
	open IN,$InFile;
	while(<IN>) {
		chomp;
		my @c=split;
		@c[3,4]=($c[3]<$c[4])? @c[3,4]:@c[4,3];
		push @{$genes_p->{"$c[1]$c[2]"}},[@c];
	}
	close IN;	
}

sub find_overlap {
	my ($c_p,$g_p,$cutoff,$tag)=@_;
	my $gene_len=abs($g_p->[4]-$g_p->[3])+1;
	my $cutoff_len;
	my $overlap=0;
	foreach my $c_gene(@$c_p) {
		my $c_gene_start=$c_gene->[3];
		my $c_gene_end=$c_gene->[4];
		my $c_gene_len=$c_gene_end-$c_gene_start+1;
		$cutoff_len=($gene_len<$c_gene_len)?$cutoff*$gene_len:$cutoff*$c_gene_len;

		if ($g_p->[4]<$c_gene_start||$g_p->[3]>$c_gene_end) {
			$overlap=0;
		}else{
			my $len=0;
			if($g_p->[4]<=$c_gene_end){
				$len=$gene_len;
			}else{
				$len=$c_gene_end-$g_p->[3]+1;
			}
			if ($len>=$cutoff_len){
				$overlap=1;
				last;
			}
		}
	}
	return $overlap;
}

sub find_nonredundance {
	my $cluster_p=shift;
	my $result_p = shift;

	my $coor = 6;

	my $max=-1;
	my $max_gene;

	foreach my $gene (@$cluster_p) {
		if ($gene->[$coor] < $max) {
			next;
		}else{
			$max=$gene->[$coor];
			$max_gene=$gene;
		}
	}
	push @{$result_p}, $max_gene;
	
	## remove the genes which overlap with the best one
	my @new_cluster; 
	my @max_temp;
	push @max_temp, $max_gene;
	foreach my $gene (@$cluster_p){
		my $overlap;
		$overlap=find_overlap(\@max_temp,$gene,$Overlap_percent,'percent');
		if ($overlap!=1) {
			push @new_cluster, $gene;
		}
	}

	&find_nonredundance(\@new_cluster, $result_p) if @new_cluster;
}
