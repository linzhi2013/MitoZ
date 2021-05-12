This is the annotation package for insect mitochondrial gene annotation. 
To run the pipeline, just type: perl MT_annotation_BGI.pl
There is no need to compile the scipt. Two perl package is needed: FindBin and File::Basename.

There is a test dataset of Drosophila melanogaster in the folder test_data for a quick run. You could follow the command line in the shell scitp: sample_command_line.sh.
perl ../MT_annotation_BGI.pl Drosophila_melanogaster_MT.fas ../MT_database/MitoProtein_RefSeq_110502.fasta ./

Then you will get a bunch of files:

# The original scaffold file and the built BLAST library:
Drosophila_melanogaster_MT.fas
Drosophila_melanogaster_MT.fas.nhr
Drosophila_melanogaster_MT.fas.nin
Drosophila_melanogaster_MT.fas.nsq
formatdb.log

# the script for tBLASTn
MitoProtein_RefSeq_110502.fasta.tblastn.shell

# the results of tBLASTn alignment
MitoProtein_RefSeq_110502.fasta.blast

# filtered results of tBLASTn
MitoProtein_RefSeq_110502.fasta.blast.filter

# the script for solar
MitoProtein_RefSeq_110502.fasta.solar.shell

# reformated tBLASTn results
MitoProtein_RefSeq_110502.fasta.blast.solar

# the filtered solar results
MitoProtein_RefSeq_110502.fasta.blast.solar.filter

# reformated solar results for clustering
MitoProtein_RefSeq_110502.fasta.blast.solar.filter.table

# non-redundant solar files
MitoProtein_RefSeq_110502.fasta.blast.solar.filter.table.nonredundance

# the corresponding information of the non-redundant hits
MitoProtein_RefSeq_110502.fasta.blast.solar.filter.nr

# the scripts for genewise annotation
MitoProtein_RefSeq_110502.fasta.genewise.shell

# the genewise results
MitoProtein_RefSeq_110502.fasta.genewise

# length of each sequences in references
MitoProtein_RefSeq_110502.fasta.length

# gff format of annotation results
MitoProtein_RefSeq_110502.fasta.solar.genewise.gff

# the corresponding nucleotide and protein sequences of annotated genes, the final annotated genes
MitoProtein_RefSeq_110502.fasta.solar.genewise.gff.cds
MitoProtein_RefSeq_110502.fasta.solar.shell


The mitochondrial database is the RefSeq MT proteins of different genes. You could update this database and just change the database file in the command line.
The folder "blast" and "solar" contains the software of BLAST and SOLAR. You could also download new versions and replace the old one. Please keep in mind that the new version of software must have the same name as old ones.
