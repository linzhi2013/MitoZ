Installation of MitoZ

Guanliang Meng, 2018-11-18

********************************************************************

# 0. System requirment: Linux

developed under: CentOS release 6.9 (Final), 2.6.32-696.30.1.el6.x86_64

# 1. Install Anaconda or Miniconda
1. Anaconda: https://anaconda.org/anaconda/python
2. Miniconda: https://conda.io/miniconda.html (recommended)


# 2. Install dependency with `conda` command

## 2.1 Set up channels

    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge

### tips
The download speeds from the above channels can be very slow, if this is the case to you,
you can set up other mirror channels, but make sure the `conda-forge` channel has the highest
priority (the last one to be added), followed by `bioconda` channel.

For example, in China, you can set:
    
    conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/
    conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/

Or,

    conda config --add channels http://mirrors.ustc.edu.cn/anaconda/cloud/bioconda/
    conda config --add channels http://mirrors.ustc.edu.cn/anaconda/cloud/conda-forge/

but I cannot ensure using mirror channels will always work. Good luck!

## 2.2 Set up an isolated enviroment for MitoZ

It is a good idea to install MitoZ into an isolated enviroment, e.g., `mitozEnv`.

    conda create  -n mitozEnv libgd=2.2.4 python=3.6.0 biopython=1.69 ete3=3.0.0b35 perl-list-moreutils perl-params-validate perl-clone circos=0.69 perl-bioperl blast=2.2.31  hmmer=3.1b2  bwa=0.7.12 samtools=1.3.1 infernal=1.1.1 tbl2asn openjdk

# 3. Activate the `mitozEnv` environment

    source activate mitozEnv

# 4. Install NCBI taxonomy database for ete3 package
1. Network connection required.
2. `HOME` directory must have more than 500M space available. If not, please refer to `https://github.com/linzhi2013/taxonomy_ranks/blob/master/README.md`

In the terminal, type `python3` then `Enter`, you will be into the Python interactive interface. Now run the following commands line-by-line in the Python interactive interface.

    from ete3 import NCBITaxa
    ncbi = NCBITaxa()
    ncbi.update_taxonomy_database()

For more details, please refer to http://etetoolkit.org/docs/latest/tutorial/tutorial_ncbitaxonomy.html


# 5. Download the MitoZ package

From `https://github.com/linzhi2013/MitoZ`

    tar -jxvf release_MitoZ_v1.0.tar.bz2
    cd release_MitoZ_v1.0
    python3 MitoZ.py

# 6. Important: make sure you are in the `mitozEnv` environment when you run MitoZ!
If you write the run commands into a script file (e.g. `work.sh`), you should also add `source activate mitozEnv` into the
script file ahead of the MitoZ commands.

********************************************************************
