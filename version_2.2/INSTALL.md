Installation of MitoZ

Guanliang Meng, 2019-01-16

********************************************************************

We provide two options to use MitoZ:
* A ready to run version of [Singularity container](https://www.sylabs.io/) [![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/1955)
* Source code with dependency to be installed with Anaconda or Miniconda.

# 1 Singularity container

## 1.1 Install Singularity
See [https://www.sylabs.io/docs/](https://www.sylabs.io/docs/) for instructions to install Singularity.

## 1.2 Download the MitoZ container

    $ singularity pull  --name MitoZ.simg shub://linzhi2013/MitoZ:v2.2

For China users, it can be difficult to pull the Singularity container for known network problem, you can
download the containers from https://pan.genomics.cn/ucdisk/s/uMFvum.

## 1.3 Run the MitoZ container

    $ /path/to/MitoZ.simg --help

## 1.4 You can also `shell` into the container

In the host, run

    $ singularity shell /path/to/MitoZ.simg

In the container, run

    $ /app/anaconda/bin/python3 /app/release_MitoZ_v2.2/MitoZ.py -h

some useful scripts are in `/app/release_MitoZ_v2.2/useful_scripts`

    $ ls -lhrt /app/release_MitoZ_v2.2/useful_scripts

To learn more about how to use Singularity, please refer to https://www.sylabs.io/docs/.

# 2 Install from source code

## 2.1 System requirment: Linux

developed under: CentOS release 6.9 (Final), 2.6.32-696.30.1.el6.x86_64

## 2.2 Install Anaconda or Miniconda
1. Anaconda: https://anaconda.org/anaconda/python
2. Miniconda: https://conda.io/miniconda.html (recommended)


## 2.3 Install dependency with `conda` command

### 2.3.1 Set up channels

    $ conda config --add channels defaults
    $ conda config --add channels bioconda
    $ conda config --add channels conda-forge

### tips
The download speeds from the above channels can be very slow, if this is the case to you,
you can set up other mirror channels, but make sure the `conda-forge` channel has the highest
priority (the last one to be added), followed by `bioconda` channel.

For example, in China, you can set:

    $ conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/
    $ conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/

Or,

    $ conda config --add channels http://mirrors.ustc.edu.cn/anaconda/cloud/bioconda/
    $ conda config --add channels http://mirrors.ustc.edu.cn/anaconda/cloud/conda-forge/

but I cannot ensure using mirror channels will always work. Good luck!

### 2.3.2 Set up an isolated enviroment for MitoZ

It is a good idea to install MitoZ into an isolated enviroment, e.g., `mitozEnv`.

    $ conda create  -n mitozEnv libgd=2.2.4 python=3.6.0 biopython=1.69 ete3=3.0.0b35 perl-list-moreutils perl-params-validate perl-clone circos=0.69 perl-bioperl blast=2.2.31  hmmer=3.1b2  bwa=0.7.12 samtools=1.3.1 infernal=1.1.1 tbl2asn openjdk

## 2.4 Activate the `mitozEnv` environment

    $ source activate mitozEnv

## 2.5 Install NCBI taxonomy database for ete3 package
1. Network connection required.
2. `HOME` directory must have more than 500M space available. If not, please refer to `https://github.com/linzhi2013/taxonomy_ranks/blob/master/README.md`

In the terminal, type `python3` then `Enter`, you will be into the Python interactive interface. Now run the following commands line-by-line in the Python interactive interface.

    from ete3 import NCBITaxa
    ncbi = NCBITaxa()
    ncbi.update_taxonomy_database()

For more details, please refer to http://etetoolkit.org/docs/latest/tutorial/tutorial_ncbitaxonomy.html


## 2.6 Download the MitoZ package

Download from `https://github.com/linzhi2013/MitoZ/tree/master/version_2.2`.

    $ tar -jxvf release_MitoZ_v2.2.tar.bz2
    $ cd release_MitoZ_v2.2
    $ python3 MitoZ.py

## 2.7 Important: make sure you are in the `mitozEnv` environment when you run MitoZ!
If you write the run commands into a script file (e.g. `work.sh`), you should also add `source activate mitozEnv` into the
script file ahead of the MitoZ commands.

********************************************************************
