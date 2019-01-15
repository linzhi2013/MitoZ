#!/usr/bin/perl -w
use strict;
use Cwd qw(abs_path);
use FindBin qw($Bin);
use Bio::SeqIO;
use Bio::SeqFeatureI;
use Getopt::Long;
use Statistics::Descriptive;

=head1 Description
	This script is an utility for chloroplast/mitochondria 
	annotation pipoline to draw circos pictures 

=head1 Author :
	YangChentao yangchentao@genomics.cn 2017/05/08

=head1 Options
	*--gb 	<str>	*.mitogenome.gb
	*--conf	<str>	*.conf
	--help		display this help information

=head1 Usage
	perl draw_circos_for_mitogenome_auto_depth.pl  -gb <*.mitogenome.gb> -conf <*.conf> 

    configure example: 
    $Bin/mitogenome.auto_depth.conf.txt

=head1 Attention
	~ If you want to draw depth part, you need set in two ways:
	1. set "depth_file = xxx" ;
	2. set "run_map = yes" and "fq = /yourpath/..."; 
	if set 1,2 in the same time, program will use 1.

=cut


my (%conf,$help,$gbfile,$configure,$warn,$outdir);

GetOptions(
			'help' => \$help,
			'gb=s' => \$gbfile,
			'conf=s' => \$configure
	);

die `pod2text $0` if ($help || !$gbfile );

$warn = <<_WARN_;
#-----------------------------------------------------------------------------------------------
WARNNING:
         No configure file!
         Here is a sample: $Bin/mitogenome.auto_depth.conf
#-----------------------------------------------------------------------------------------------
_WARN_


die  "$warn" unless ($configure);

# reading configures from $configure
open CC, "$configure";	# circos_path,win,cds,rRNA,tRNA,locus_color,label_color,gc_fill,depth_fill
while (<CC>){
	next if (/^#/);
	next if (/^\s*$/);
	chomp;
	my @cc = split /=/,$_;
	$cc[0] =~ s/\s//g;
	if (/opts_samtools/) {
		$cc[1] =~ s/^\s+//;
		$cc[1] =~ s/\s+$//;
	}else{
		$cc[1] =~ s/\s//g;
		s/\s//g;
	}
	$conf{$cc[0]} = $cc[1];

}

close CC;
$outdir = $conf{'outdir'};
$outdir = abs_path($outdir);

## check configure
die "circos path is bad !" if ( ! $conf{'circos_path'} ) ;
# outdir 
`[ -d $outdir ] || mkdir $outdir`;

$conf{'win'} = 50 if (! $conf{'win'}); # slipping windows size for GC content, default=50 bp
$conf{'gc'} = 'yes' if (! $conf{'gc'});
$conf{'run_map'} = 'yes' if (! $conf{'run_map'});
$conf{'base'} = 'no' if (! $conf{'base'});
$conf{'threads'} ||= 2;

$conf{'opts_samtools'} ||= '' ; # optional arguments for samtools.
my $opts_samtools;
if ($conf{'opts_samtools'}){
	$opts_samtools = $conf{'opts_samtools'};

}

# calculate depth
my $bwa = $conf{'bwa'} ;
my $samtools = $conf{'samtools'};
# turn on mapping function
my $mapping_on = 0;
my $generated_depth ;


if ($conf{'depth_file'} and $conf{'run_map'} eq 'yes'){
	
	print "WARNNING: you set both depth_file and run_map = yes, so I will just use existed depth file\n";
	
	$generated_depth = $conf{'depth_file'};


}elsif ($conf{'run_map'} eq 'yes') {
		
	$mapping_on = 1;
	die "bwa path is bad!" if ( ! $bwa ) ;
	die "samtools path is bad!" if ( ! $samtools ) ;
	if( ! $conf{'fq'}) {
		die "Fastq file is necessary in configures when you set \"run_map = yes\"" ;
	}else{
		open FA,">$outdir/circos.mito.fa";
		open MAP, ">$outdir/circos.map.sh";
		my @fq = split /,/,$conf{'fq'} ;
		my ($fq1, $fq2);
		if (@fq < 2){
			die "fastq_1 is bad file!" unless (-e $fq[0]) ;
			$fq1 = abs_path($fq[0]);
			$fq2 = ""
		}else{
			die "fastq_1 is bad file!" unless (-e $fq[0]) ;
			die "fastq_2 is bad file!" unless (-e $fq[1]) ;
			$fq1 = abs_path($fq[0]);
			$fq2 = abs_path($fq[1]);
		}

		print MAP "$bwa index $outdir/circos.mito.fa \n";
		print MAP "$bwa mem -t $conf{'threads'} $outdir/circos.mito.fa $fq1  $fq2 |samtools view -bS -q 30 -h -o $outdir/circos.bam\n";
		print MAP "$samtools sort $outdir/circos.bam -o $outdir/circos.sorted.bam\n";
		print MAP "$samtools depth $opts_samtools $outdir/circos.sorted.bam > $outdir/circos.dep\n";
		print MAP "awk \'{print \$1\,\$2,\$2,\$3}\' $outdir/circos.dep > $outdir/circos.depth.txt\n";

	}

}elsif ($conf{'depth_file'}) {

	$generated_depth = $conf{'depth_file'} ;
}


if ($conf{'base'} eq 'yes') {
	open BASE,">$outdir/circos.base.txt";
}

if ($conf{'gc'} eq 'yes') {
	open GC,">$outdir/circos.fa.gc.txt";
}

### read genebank file 
# whether link chromsome as a whole according to topology[linear|circular] of mitogenome.
my %topology;
my $locus_line = `grep 'LOCUS' $gbfile `;
chomp $locus_line;
my @lines = split /\n/,$locus_line;
my $submitos = @lines;
if ($submitos > 1){
	print "WARNNING: $submitos sequences will be presented on circos figure.\n";
}


# open file handles
open KAR,">$outdir/circos.karyotype.txt";
open FEA,">$outdir/circos.features.txt"; # location information
open TEXT,">$outdir/circos.gene.text.txt"; # gene name text

open ORI,">$outdir/circos.strand.text.txt"; # strand label
print ORI "mt1\t0\t300\t+\tr0=1r-150p,r1=1r-100p\n";
close ORI;

my $break;
foreach my $l(@lines) {
	my $a = (split/\s+/,$l)[1];
	my $topo = (split/\s+/,$l)[5];
	$topology{$a} = $topo;
	if ($topo eq 'circular') {
		$break = 0;
	}else {
		$break = "0.5r";
	}
}
# get cds rRNA tRNA's location and strand information
my $in = Bio::SeqIO-> new(-file => "$gbfile", "-format" => 'genbank');

my $id = 0;
my ($chr,$mt);
my %locus_id;

while (my $seq_obj=$in->next_seq()) {

	my $source = $seq_obj->seq;
	my $sou_len = $seq_obj->length;
	my $locus = $seq_obj -> display_id;

	$id ++;
	$chr = "chr" . $id;
	$mt = "mt" . $id;
	$locus_id{$locus} = $mt;
	
	#print KAR "$chr - $mt\t$locus-$topology{$locus}\t0\t$sou_len\tgrey\n";
	print KAR "$chr - $mt\t$locus\t0\t$sou_len\tgrey\n";
	
	print FA ">$mt\n$source\n" if ($mapping_on) ;

	if ($conf{'base'} eq 'yes') {
		my @base = split//,$source;
		# output base 
		foreach my $i(0..$#base) {
			my $j = $i +1;
			print BASE "$mt\t$j\t$j\t$base[$i]\n";
		}
	}
	
	#calculate GC content
	if ($conf{'gc'} eq 'yes') {
		my $win = $conf{'win'};
		for (my $i = 0; $i < $sou_len - $win -1; $i+=$win) {
			my $tmp = substr($source,$i,$win);
			my $gc_num = $tmp =~ tr/GCgc/GCgc/;
			my $GC = $gc_num/$win;
			my $start = $i +1 ;
			my $end = $i + $win ;
			print GC "$mt\t$start\t$end\t$GC\n";
		}	
	}

	#print ORI "$mt\t0\t300\t+\tr0=1r-200p,r1=1r-100p\n";
	#print ORI "$mt\t0\t300\t-\tr0=1r+100p,r1=1r+150p\n";
	
	for my $feature ($seq_obj->top_SeqFeatures){
		my ($db_xref,$val,$location);

		if ($feature->primary_tag eq 'CDS' )  {
			my $seq = $feature->spliced_seq->seq;
			my $start = $feature -> start;
			my $end = $feature -> end;
			my $strand = $feature -> strand;
			my $direction;
			if ($strand == 1) {
				$direction = '+';
				print FEA "$mt\t$start\t$start\tfill_color=black,r0=0.965r,r1=1r\n";
				print FEA "$mt\t$end\t$end\tfill_color=black,r0=0.965r,r1=1r\n";
				print FEA "$mt\t$start\t$end\tfill_color=$conf{'cds'},r0=0.965r,r1=1r\n";
			}elsif ($strand == -1){
				$direction = '-';
				print FEA "$mt\t$start\t$start\tfill_color=black,r0=1r,r1=1.035r\n";
				print FEA "$mt\t$end\t$end\tfill_color=black,r0=1r,r1=1.035r\n";
				print FEA "$mt\t$start\t$end\tfill_color=$conf{'cds'},r0=1r,r1=1.035r\n";
			}else {
				$direction = '?';
				print FEA "$mt\t$start\t$start\tfill_color=black,r0=0.9825r,r1=1.0175r\n";
				print FEA "$mt\t$end\t$end\tfill_color=black,r0=0.9825r,r1=1.0175r\n";
				print FEA "$mt\t$start\t$end\tfill_color=$conf{'cds'},r0=0.9825r,r1=1.0175r\n";

			}
			#print CDS "$mt\t$start\t$end\n";

			if ($feature->has_tag('gene')) {
				for $val ($feature->get_tag_values('gene')){
				
					print TEXT "$mt\t$start\t$end\t$val\n";
				}
			}else{
			
				print TEXT "$mt\t$start\t$end\tCDS_NA\n";
			}
		}elsif ($feature->primary_tag eq 'rRNA' )  {
		
			my $start = $feature -> start;
			my $end = $feature -> end;
		#	print RRNA "$mt\t$start\t$end\n";

			my $strand = $feature -> strand;
			my $direction;
			if ($strand == 1) {
				$direction = '+';
				print FEA "$mt\t$start\t$start\tfill_color=black,r0=0.965r,r1=1r\n";
				print FEA "$mt\t$end\t$end\tfill_color=black,r0=0.965r,r1=1r\n";
				print FEA "$mt\t$start\t$end\tfill_color=$conf{'rRNA'},r0=0.965r,r1=1r\n";
			}elsif ($strand == -1){
				$direction = '-';
				print FEA "$mt\t$start\t$start\tfill_color=black,r0=1r,r1=1.035r\n";
				print FEA "$mt\t$end\t$end\tfill_color=black,r0=1r,r1=1.035r\n";
				print FEA "$mt\t$start\t$end\tfill_color=$conf{'rRNA'},r0=1r,r1=1.035r\n";
			}else {
				$direction = '?';
				print FEA "$mt\t$start\t$start\tfill_color=black,r0=0.9825r,r1=1.0175r\n";
				print FEA "$mt\t$end\t$end\tfill_color=black,r0=0.9825r,r1=1.0175r\n";
				print FEA "$mt\t$start\t$end\tfill_color=$conf{'rRNA'},r0=0.9825r,r1=1.0175r\n";
			}
			if ($feature->has_tag('gene')) {
				for $val ($feature->get_tag_values('gene')){
					print TEXT "$mt\t$start\t$end\t$val\n";
				}
			}else{
			
				print TEXT "$mt\t$start\t$end\trRNA_NA\n";
			}

		}elsif ($feature->primary_tag eq 'tRNA' )  {
		 
		 	my $start = $feature -> start;
			my $end = $feature -> end;
		#	print TRNA "$mt\t$start\t$end\n";

			my $strand = $feature -> strand;
			my $direction;
			if ($strand == 1) {
				$direction = '+';
				print FEA "$mt\t$start\t$start\tfill_color=black,r0=0.965r,r1=1r\n";
				print FEA "$mt\t$end\t$end\tfill_color=black,r0=0.965r,r1=1r\n";
				print FEA "$mt\t$start\t$end\tfill_color=$conf{'tRNA'},r0=0.965r,r1=1r\n";
			}elsif ($strand == -1){
				$direction = '-';
				print FEA "$mt\t$start\t$start\tfill_color=black,r0=1r,r1=1.035r\n";
				print FEA "$mt\t$end\t$end\tfill_color=black,r0=1r,r1=1.035r\n";
				print FEA "$mt\t$start\t$end\tfill_color=$conf{'tRNA'},r0=1r,r1=1.035r\n";
			}else {
				$direction = '?';
				print FEA "$mt\t$start\t$start\tfill_color=black,r0=0.9825r,r1=1.0175r\n";
				print FEA "$mt\t$end\t$end\tfill_color=black,r0=0.9825r,r1=1.0175r\n";
				print FEA "$mt\t$start\t$end\tfill_color=$conf{'tRNA'},r0=0.9825r,r1=1.0175r\n";
			}
			if ($feature->has_tag('gene')) {
				for $val ($feature->get_tag_values('gene')){
					
					print TEXT "$mt\t$start\t$end\t$val\n";
				}
			}else{	
					
					print TEXT "$mt\t$start\t$end\ttRNA_NA\n";
			}
		}
	}
}

# close file handles

if ($mapping_on && $conf{'fq'}) {
	close FA;
	close MAP;
	print STDERR "run mapping shell\n";
	system("sh $outdir/circos.map.sh");
	`rm $outdir/circos.mito.fa.*`;
	unless(-s "$outdir/circos.depth.txt"){
		print "Something wrong with mapping process!\n";
		exit();
	}

	$generated_depth = "$outdir/circos.depth.txt";

}elsif ($conf{'depth_file'}){
	$generated_depth = &format_depth_file($generated_depth);
}

if ($conf{'base'} eq 'yes') {
	close BASE;
}
if ($conf{'gc'} eq 'yes') {
	close GC;
}

close KAR;
close FEA;
close TEXT;


# creat a circos.conf file
open CON,">$outdir/circos.conf";
print  CON <<_CONFIG_;

<<include etc/colors_fonts_patterns.conf>>

#-----------------image------------------
<image>
###<<include etc/image.conf>>
dir   = $conf{'outdir'}
file  = circos.png
png   = $conf{'png'}
svg   = $conf{'svg'}

# radius of inscribed circle in image
radius         = 1500p

# by default angle=0 is at 3 o'clock position
angle_offset      = -90

#angle_orientation = counterclockwise
auto_alpha_colors = yes
auto_alpha_steps  = 5
background = $conf{'background'}
</image>

#-----------------ideogram------------------
<ideogram>

<spacing>
default = 0.01r
break   = $break
</spacing>

###<<include ideogram.position.conf>>
radius           = 0.82r
thickness        = 20p
fill             = yes
fill_color       = grey
stroke_thickness = 3
stroke_color     = black

###<<include ideogram.label.conf>>
show_label       = yes
label_font       = bolditalic
label_radius     = dims(ideogram,radius_outer) - 0.1r
#label_radius     = 0.2r
label_size       = 28
label_parallel   = yes
label_case       = lower
#label_format     = eval(sprintf("chr%s",var(label)))
#label_format     = eval(var(labe))


###<<include bands.conf>>
show_bands            = yes
fill_bands            = yes
band_stroke_thickness = 2
band_stroke_color     = white
band_transparency     = 0

</ideogram>
#-----------------ticks------------------
show_ticks          = yes
show_tick_labels    = yes

<ticks>

radius           = dims(ideogram,radius_outer)
#radius           = 1r+0.06r
orientation      = out
label_multiplier = 1e-3
color            = black
thickness        = 2p
font             = blod

<tick>
spacing        = 1u
show_label     = yes
label_size     = 25p
size           = 25p
format         = %d
label_offset   = 2p
#suffix         = " kb"
</tick>

<tick>
spacing        = 5u
show_label     = yes
label_size     = 30p
size           = 30p
format         = %d
suffix         = " kb"
label_offset   = 2p
</tick>

<tick>
spacing        = 10u
show_label     = yes
label_size     = 30p
size           = 30p
format         = %d
label_offset   = 2p
suffix         = " kb"
</tick>

</ticks>
#-----------------karyotype------------------

karyotype   = circos.karyotype.txt

chromosomes_units = 1000
#chromosomes       = mt1
chromosomes_display_default = yes

#-----------------plots------------------

<plots>

############ gene name text
<plot>

type       = text
color      = $conf{'label_color'}
label_font = default
label_size = 28p
file = circos.gene.text.txt
r1   = 1r+300p
r0   = 1r+10p
show_links     = yes
link_dims      = 0p,0p,70p,0p,10p
link_thickness = 2p
link_color     = red

label_snuggle         = yes
max_snuggle_distance  = 1r
snuggle_tolerance     = 0.25r
snuggle_sampling      = 2

</plot>

<plot>

type       = text
color      = $conf{'label_color'}
label_font = bold
label_size = 40p
file = circos.strand.text.txt
show_links     = no

</plot>


_CONFIG_

if ($conf{'gc'} eq 'yes') {
	print CON <<_CONFIG_;

###############GC content
<plot>
type      = histogram
file      = circos.fa.gc.txt

r1        = 0.615r
r0        = 0.45r
max       = 1
min       = 0

stroke_type = line
thickness   = 2
color       = $conf{'gc_fill'}
extend_bin  = no
fill_color = $conf{'gc_fill'}
<backgrounds>
#<background>
#y1    = -0.1
#color = lred
#</background>
#<background>
#y0    = 0
#color = lgreen
#</background>
</backgrounds>

<axes>

<axis>
spacing   = 0.05r
color     = lgrey
thickness = 1
</axis>

<axis>
position = 0.5r
color = dred
thickness = 2
</axis>  

</axes>

<rules>
use = no
<rule>
condition  = var(value) >0.50
fill_color = dyellow
</rule>

</rules>

</plot>

_CONFIG_

}

if ($conf{'base'} eq 'yes') {
	print CON <<_CONFIG_;

########## sequence base
<plot>
type       = text
label_font = mono
file       = circos.base.txt
r1         = 0.91r
r0         = 0.88r
label_size = 20
padding    = -0.25r
#label_rotate = no

<rules>
<rule>
condition = var(value) eq "A"
color     = red
</rule>
<rule>
condition = var(value) eq "T"
color     = blue
</rule>
<rule>
condition = var(value) eq "C"
color     = green
</rule>
</rules>
</plot>

_CONFIG_
}



if ($generated_depth) {

	# caculate upper quartile value
	my $str = `grep -v "^#" $generated_depth|awk '{print \$4}'`;
	chomp $str;
	my @tmp = split /\n/,$str;
	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@tmp);	
	my $upper_quartile = $stat->quantile(3);
	my $max_depth = $stat->quantile(4);

	print CON <<_CONFIG_;

########### depth information
<plot>
type      = line
thickness = 2
max_gap = 1u
file    = $generated_depth
color   = dgreen
min     = 0
max     = $max_depth
r0      = 0.618r
r1      = 0.768r
fill_color = $conf{'depth_fill'}


<axes>
<axis>
color     = lgrey_a2
thickness = 1
spacing   = 0.06r
</axis>
</axes>

<rules>

<rule>
condition    = var(value) > $upper_quartile
color        = $conf{'depth_fill'}
fill_color   = $conf{'depth_fill'}
</rule>

<rule>
condition    = var(value) < 20
color        = dred
fill_color   = dred_a1
</rule>

</rules>
</plot>

_CONFIG_

}

print CON "</plots>\n";

print CON <<_CONFIG_;

#-----------------highlights------------------
<highlights>

# CDS & rRNA & tRNA
<highlight>
file         = circos.features.txt
</highlight>

</highlights>

<<include etc/housekeeping.conf>>

_CONFIG_

close CON;


`$conf{'circos_path'} -conf $outdir/circos.conf`;
print "All done!\n";
print "PNG => $outdir/circos.png\n";
print "SVG => $outdir/circos.svg\n";


sub format_depth_file{
	open IN,shift;
	my $output = "$outdir/tmp.depth.txt";
	open OUT, ">$output";
	while (<IN>) {
		next if (/^#/);
		chomp;
		my @tmp = split /\s+/;
		if (exists $locus_id{$tmp[0]}) {
			print OUT "$locus_id{$tmp[0]}\t$tmp[1]\t$tmp[1]\t$tmp[2]\n";
		}else{
			print "bad depth_file: $tmp[0] does not in genebank file!\n";
			exit();
		}

	}
	close IN;
	close OUT;
	`mv $output $outdir/circos.depth.txt`;
	return "$outdir/circos.depth.txt";
}

