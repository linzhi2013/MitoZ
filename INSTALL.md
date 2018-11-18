Installation of MitoZ

Guanliang Meng, 2018-11-18

********************************************************************

# 0. System requirment: Linux

developed under: CentOS release 6.9 (Final), 2.6.32-696.30.1.el6.x86_64

# 1. Install Anaconda or Miniconda
1. Anaconda: https://anaconda.org/anaconda/python
2. Miniconda: https://conda.io/miniconda.html


# 2. Install dependency with `conda`

## 2.1 Set up channels

    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge

## 2.2 Set up an isolated enviroment for MitoZ

It is a good idea to install MitoZ into an isolated enviroment, e.g., `mitozEnv`.

    conda create  -n mitozEnv  -c conda-forge libgd=2.2.4 python=3.6.0 biopython==1.69 ete3==3.0.0b35
    conda install -y -n mitozEnv -c bioconda perl-list-moreutils perl-params-validate perl-clone circos=0.69 perl-bioperl


# 3. Activate the `mitozEnv` environment

    source activate mitozEnv

# 4. Install NCBI taxonomy database for ete3 package
1. Network connection required.
2. `HOME` directory must have more than 500M space available.

In the terminal, type `python3` then `Enter`, you will be into the Python interactive interface. Now run the following commands line-by-line in the Python interactive interface.

    from ete3 import NCBITaxa
    ncbi = NCBITaxa()
    ncbi.update_taxonomy_database()

For more details, please refer to http://etetoolkit.org/docs/latest/tutorial/tutorial_ncbitaxonomy.html

********************************************************************
