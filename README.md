# Manual of MitoZ

# 1. About MitoZ

MitoZ is a Python3-based toolkit which aims to automatically filter pair-end raw data (fastq files), assemble genome, search for mitogenome sequences from the genome assembly result, annotate mitogenome (genbank file as result), and mitogenome visualization. MitoZ is available from `https://github.com/linzhi2013/MitoZ`.


# 2. System requirment

## 2.1 Platform

MitoZ is developed and tested under `Linux version 2.6.32-696.el6.x86_64 (mockbuild@c1bm.rdu2.centos.org) (gcc version 4.4.7 20120313 (Red Hat 4.4.7-18) (GCC) ) #1 SMP Tue Mar 21 19:29:05 UTC 2017`.

## 2.2 Harddisk space

\>100G

## 2.3 Memory

It takes ~100G when we tested MitoZ with `--thread_number 16`. This is because MitoZ uses the *de Bruijn* graph (DBG) algorithm to perform *de novo* assembly. To learn more, see https://doi.org/10.1093/bioinformatics/btu077.


# 3. Get started

MitoZ includes multiple modules, including `all`, `all2`, `filter`, `assemble`, `findmitoscaf`, `annotate` and `visualize`.

**Important: make sure you are in the `mitozEnv` environment when you run MitoZ!**

    source activate mitozEnv

see `INSTALL.md` for installation instruction.


Now create a directory for one sample:

    mkdir ~/example
    cd ~/example

# 4. Data requirement

## 4.1 Material source

The most preferable data for mitochondrial genome assembly is always data with high
ratio of mitochondrial derived reads and little contamination. For example, tissue samples
may be better than blood samples or gut samples.

## 4.2 Requirment of data size and insert size

About 1.5 to 3G base pair (bp) is enough for mitochondrial genome assembly.

e.g. `raw.1.fq.gz` and `raw.2.fq.gz` or `clean.1.fq.gz` and `clean.2.fq.gz`

The read length should be >= 71bp (PE71). Typically, I use data of PE100 or PE150 sequencing.

The length of read1 and read2 must be equal. You should trim your data (e.g. use the option `--keep_region` in `filter` function of MitoZ) before running MitoZ.

The insert size of pair-end library should be small insert size (<1000 bp). I do not recommend to use data of large insert size (>1000bp, mate-pair library), because this kind of data usually is not good as small insert size data.


## 4.3 Data pretreatment

MitoZ supports simple data pretreatment (remove low quality, many Ns reads, duplications),
thus you can provide MitoZ the raw data (fastq files) from WGS experiments directly.
In this case, you can use the `all` or `filter` module to perform data filtering.

Or, you can provide MitoZ the clean data, which have been filtered by other tools.

## 4.4 Modules support both single-end data and pair-end data

* `all2`

* `assemble`

* `annotate`

* `visualize`


## 4.5 Modules support pair-end data only

* `all`

* `filter`

* `findmitoscaf`(needs fastq only when the input assembly (containing nuclear and mitochondrial sequences) is not from SOAPdenovo-Trans or mitoAssemble, to caculate the sequence sequencing depth)


## 4.6 Fasta file

When you annotate a mitogenome sequence(s) stored in fasta file, the sequence id can not be too long, or MitoZ will fail. This is for some limitation in BioPython that MitoZ invokes.


## 4.7 Genetic code

It is important to set a correct genetic code for MitoZ (`--genetic_code` option). Usually, arthropods use the invertebrate mitochondrial code (`--genetic_code 5`), and mammals use the vertebrate mitochondrial code (`--genetic_code 2`).Please refer to https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for more details.

MitoZ use the genetic code for annotating the Protein coding genes (PCGs) of mitochondrial genome.


## 4.8 Taxa group

It is also important to set a correct taxa group for MitoZ (`--clade` option).

MitoZ use the `--clade` option to choose corresponding database (HMM modules, CM modules, protein reference sequences).

# 5. Supporting for specifying parameters in a configure file

Now MitoZ supports use a configure file to set the parameter, besides using the command
line way (as the examples in the later sections).

## 5.1 create a configure file for desired module

For example,

    python3 MitoZ.py all --create_config

will a file `mitoz_all_config.txt`, which contains the paramters same as output by `python3 MitoZ.py all `.

Then modify the file `mitoz_all_config.txt` as the instructions in the file.


## 5.2 run the module

    python3 MitoZ.py all --config mitoz_all_config.txt


# 6. `all` module

`all` module supports pair-end data only.

`all` module requires only two input pair-end fastq files, and outputs a genbank file containing mitochondrial genome sequences and annotation information.

## 6.1 Input files

Pair-end(PE) fastq files (`raw.1.fq.gz` and `raw.2.fq.gz`), and optional files `1.adapter.list.gz` and `2.adapter.list.gz`.

## 6.2 Example

    python3 MitoZ.py all --genetic_code 5 --clade Arthropoda --outprefix ZZZ \
    --thread_number 12 \
    --fastq1 raw.1.fq.gz \
    --fastq2 raw.2.fq.gz \
    --fastq_read_length 150 \
    --insert_size 250  \
    --run_mode 2 \
    --filter_taxa_method 1 \
    --requiring_taxa 'Arthropoda'

For more details, please refer to `python3 MitoZ.py all -h`

## 6.3 directory structure


    example
    ├── tmp
    │   ├── ZZZ.annotation
    │   ├── ZZZ.assembly
    │   └── ZZZ.cleandata
    └── ZZZ.result
        ├── work71.hmmtblout.besthit.sim.filtered.high_abundance_10.0X.reformat.sorted.Not-picked
        ├── work71.hmmtblout.besthit.sim.filtered.high_abundance_10.0X.reformat.sorted.Not-picked.fa
        ├── work71.hmmtblout.besthit.sim.filtered.low_abundance
        ├── work71.hmmtblout.besthit.sim.filtered.low_abundance.fasta
        ├── work71.mitogenome.fa
        ├── work71.most_related_species.txt
        ├── ZZZ_mitoscaf.fa.gbf
        ├── ZZZ_mitoscaf.fa.sqn
        ├── ZZZ_mitoscaf.fa.tbl
        ├── ZZZ_mitoscaf.fa.val
        ├── errorsummary.val
        ├── ZZZ.most_related_species.txt
        ├── ZZZ.fasta
        ├── ZZZ.cds
        ├── ZZZ.rrna
        ├── ZZZ.trna
        ├── circos.png
        ├── circos.svg
        ├── summary.txt
        └── README.txt

The intermediate files are in the `tmp` directory, and `ZZZ.result` contains the result files
for you sample.

Files in `ZZZ.result` including:

* `README.txt`

A brief instructions about each file in the `ZZZ.result` directory.

* `summary.txt`

A summary about the mitogenome in `ZZZ_mitoscaf.fa.gbf`, including a list of genes,
numbers of genes recovered totally, and what genes may be missing.

* `circos.png` and `circos.svg`

Visualizations of the mitogenome in `ZZZ_mitoscaf.fa.gbf`.

* `ZZZ.most_related_species.txt`

The most closely related species of the mitogenome in `ZZZ_mitoscaf.fa.gbf`. This information
comes protein coding gene (PCG) annotations, during which each PCG will be assigned to a protein
sequence of related species by the Genewise program. The most occurrence of related speices
for each sequence will be reported in `ZZZ.most_related_species.txt`.

* `ZZZ_mitoscaf.fa.gbf`, `ZZZ.fasta`, `ZZZ_mitoscaf.fa.sqn` and `ZZZ_mitoscaf.fa.tbl`.

The mitogenome files in different format. The sequences whose sequence ids have `_FivePCGs` suffixs
(if any) are not considered as our mitogenome of target species, and they are output intendedly
for further inspection by users.

* `ZZZ.cds`, `ZZZ.rrna`, `ZZZ.trna`

The individual gene sequences in fasta format.

* `ZZZ_mitoscaf.fa.val` and `errorsummary.val`

The two files are generated by NCBI's tbl2asn program, describing related warnings and errors
for file `ZZZ_mitoscaf.fa.gbf`. Thus, you should check this file to confirm if there are any
assembly or annotation problems.

For example, a InternalStop for protein coding genes can be caused by assemlby (e.g., there are Ns),
annotation, mutation, incorrect genetic code, or the sequence simply comes from nuclear
mitochondrial DNA segments (NUMTs). You may inspect sequence sequencing depth from the visualization
result files (`circos.png` and `circos.svg`).


* `work71.mitogenome.fa` and `work71.most_related_species.txt`

The mitogenome sequences from each kmer (e.g. kmer 71) assembly and its most closely related species.
Sequences in `work71.mitogenome.fa` can be different from those in `ZZZ_mitoscaf.fa.gbf` and `ZZZ.fasta`, if the mitogenome sequence in `work71.mitogenome.fa` has overlapping regions between
5' and 3' ends. In such a case, MitoZ will detect the overlapping region and trim it, resulting
files `ZZZ_mitoscaf.fa.gbf` and `ZZZ.fasta`.

* `work71.hmmtblout.besthit.sim.filtered.high_abundance_10.0X.reformat.sorted.Not-picked`, `work71.hmmtblout.besthit.sim.filtered.high_abundance_10.0X.reformat.sorted.Not-picked.fa`, `work71.hmmtblout.besthit.sim.filtered.low_abundance` and `work71.hmmtblout.besthit.sim.filtered.low_abundance.fasta`

These `*.low_abundance*` and `*.high_abundance*` files, which are the sequences with low abundances
or high abundances but not selected as outputs by MitoZ.

These sequences, together with the sequences whose sequence ids have `_FivePCGs` suffixs in `ZZZ_mitoscaf.fa.gbf` and `ZZZ.fasta` (if any), may be useful if, for example, you want to find
potential NUMTs, or want to know if there some sequences of non-target-species.



When you use other modules, some of those files (directories) may be absent.


# 7 `all2` module

`all2` module supports single-end data and pair-end data.

`all2` is amolst the same as `all`, except that `all2` doesn't filter the input fastq files.

## 7.1 Input files
Pair-end(PE) fastq files (`clean.1.fq.gz` and `clean.2.fq.gz`).

## 7.2 Example

For pair-end data:

    python3 MitoZ.py all --genetic_code 5 --clade Arthropoda --outprefix test \
    --thread_number 8 \
    --fastq1 clean.1.fq.gz \
    --fastq2 clean.2.fq.gz \
    --fastq_read_length 150 \
    --insert_size 250 \
    --run_mode 2 \
    --filter_taxa_method 1 \
    --requiring_taxa 'Arthropoda'

For single-end data:

    python3 MitoZ.py all --genetic_code 5 --clade Arthropoda --outprefix test \
    --thread_number 8 \
    --fastq1 clean.1.fq.gz \
    --fastq_read_length 150 \
    --insert_size 250 \
    --run_mode 2 \
    --filter_taxa_method 1 \
    --requiring_taxa 'Arthropoda'

# 8. `filter` module

`filter` module supports pair-end data only.

`filter` is to filter input raw fastq files (`raw.1.fq.gz` and `raw.2.fq.gz`), outputs clean fastq files (`clean.1.fq.gz` and `clean.2.fq.gz`)

## 8.1 Input files
Pair-end(PE) fastq files (`raw.1.fq.gz` and `raw.2.fq.gz`), and optional files `1.adapter.list.gz` and `2.adapter.list.gz`.

## 8.2 Example

    python3 MitoZ.py filter \
    --fastq1 raw.1.fq.gz \
    --fastq2 raw.2.fq.gz \
    --fastq3 clean.1.fq.gz \
    --fastq4 clean.2.fq.gz \
    --outprefix test


# 9. `assemble` module

`assemble` module support single-end data and pair-end data

`assemble` is to assemble `clean.1.fq.gz` and `clean.2.fq.gz`, search for mitochondrial sequences from the assembly. Output is a mitochondrial sequence file in fasta format.

## 9.1 Input files
Pair-end(PE) fastq files (`clean.1.fq.gz` and `clean.2.fq.gz`).

## 9.2 Example

For pair-end data:

    python3 MitoZ.py assemble --genetic_code 5 --clade Arthropoda --outprefix test \
    --thread_number 8 \
    --fastq1 clean.1.fq.gz \
    --fastq2 clean.2.fq.gz \
    --fastq_read_length 150 \
    --insert_size 250 \
    --run_mode 2 \
    --filter_taxa_method 1 \
    --requiring_taxa 'Arthropoda'

For single-end data:

    python3 MitoZ.py assemble --genetic_code 5 --clade Arthropoda --outprefix test \
    --thread_number 8 \
    --fastq2 clean.2.fq.gz \
    --fastq_read_length 150 \
    --insert_size 250 \
    --run_mode 2 \
    --filter_taxa_method 1 \
    --requiring_taxa 'Arthropoda'


# 10. `findmitoscaf` module

`findmitoscaf` module supports pair-end data only.

`findmitoscaf` is to search for the mitochondrial sequences from fasta file which contains non-mitochondrial sequences. Output is a mitochondrial sequence file in fasta format.

## 10.1 Input files

* `work71.scafSeq`

A file contains non-mitochondrial sequences.

If `work71.scafSeq` is generated by SOAPdenovo-Trans-127mer or mitoAssemble, you can specify the option `--from_soaptrans`, and in this case, `work71.scafSeq` is the only input file you need.

Otherwise, you still need following two files as input,

* `clean.1.fq.gz` and `clean.2.fq.gz`


## 10.2 Example

    python3 MitoZ.py findmitoscaf --genetic_code 5 --clade Arthropoda --outprefix test \
    --thread_number 8 \
    --from_soaptrans \
    --fastafile work71.scafSeq

Or,

    python3 MitoZ.py findmitoscaf --genetic_code 5 --clade Arthropoda \
    --outprefix test --thread_number 8 \
    --fastq1 clean.1.fq.gz \
    --fastq2 clean.2.fq.gz \
    --fastq_read_length 150 \
    --fastafile work71.scafSeq


# 11. `annotate` module

`annotate` module supports single-end and pair-end data.

`annotate` is to annotate the input mitogenome sequence, including protein coding genes (PCGs), tRNA genes and rRNA genes. Output is a genbank file containing mitochondrial genome sequences and annotation information.

## 11.1 Input files
e.g. `mitogenome.fa`

A fasta file containing the mitochondrial seqeunces.

## 11.2 Example

    python3 MitoZ.py annotate --genetic_code 5 --clade Arthropoda \
    --outprefix test --thread_number 8 \
    --fastafile mitogenome.fa


If you want to see the abundance along the mitogenome sequences, you will also need to set `--fastq1` and/or `--fastq2`.

# 12. `visualize` module

`visualize` module supports single-end and pair-end data.

`visualize` module is to visualize the genbank file.

## 12.1 Input files
e.g. `mitogenome.gb`

A Genbank file.

## 12.2 Example

    python3 MitoZ.py visualize --gb mitogenome.gb

If you want to show sequencing depth along the mitogenome, you need to set `--fastq1` and/or `--fastq2` and `--depth`.


# 13. Multi-Kmer mode
when there are missing PCGs after you run MitoZ in quick mode (`--run_mode 2`), you can try with the multi-Kmer mode (`--run_mode 3`).

You should provide the quick mode assembly as input, including files:

1. `work71.hmmout.fa`, or a file (e.g. `quickMode.fa`) which you convince containing the correct mitogenome sequences for your sample. And manually create a file (e.g. `quick_mode_fa_genes.txt`) describing what PCG genes on each sequence, format:

        seqid1 PCG1 PCG2
        seqid2 PCG3

2. `work71.hmmtblout.besthit.sim.filtered.fa`

3. `work71.hmmtblout.besthit.sim.filtered.high_abundance_*X.reformat.sorted`

in the directory of `outprefix.assembly`.

## 13.1 Example

    python3 MitoZ.py all2 --genetic_code 5 --clade Arthropoda --outprefix test \
    --thread_number 12 --fastq1 clean.1.fq.gz --fastq2 clean.2.fq.gz \
    --fastq_read_length 150 --insert_size 250 \
    --run_mode 3 \
    --filter_taxa_method 1 \
    --requiring_taxa 'Arthropoda' \
    --quick_mode_seq_file quickMode.fa \
    --quick_mode_fa_genes_file quick_mode_fa_genes.txt  \
    --missing_PCGs ND4L ND6 ND2 \
    --quick_mode_score_file work71.hmmtblout.besthit.sim.filtered.high_abundance_10.0X.reformat.sorted  \
    --quick_mode_prior_seq_file work71.hmmtblout.besthit.sim.filtered.fa

The result file is `outprefix.multiKmer_seq_picked.clean.fa` under directory `outprefix.assembly2`.


# 14. Rearrangement of mitogenome basing on reference mitogeome


# Useful scripts

    python3 bin/annotate/genbank_file_tool.py
    usage: genbank_file_tool.py [-h] {cut,comrev,sort,select} ...

    Description
        A tool to deal with genbank records.

    Version
        0.0.1

    Author
        mengguanliang(at) genomics (dot) cn, BGI-Shenzhen.

    positional arguments:
      {cut,comrev,sort,select}
        cut                 cutting sequences (5' and/or 3' end).
        comrev              get complement reverse of genbank records
        sort                sort the gene orders (input should all be circular
                            records!!!)
        select              output specific genbank records

    optional arguments:
      -h, --help            show this help message and exit

# change log
v2.2
201801
1. add configure file support
2. add single-end data support
3. remove the useless codes
4. to be more user-friendly


# Citation

    Guanliang Meng, Yiyuan Li, Chentao Yang, Shanlin Liu. MitoZ: A toolkit for mitochondrial genome assembly, annotation and visualization; doi: https://doi.org/10.1101/489955


