FROM ubuntu:16.04

MAINTAINER Guanliang Meng <linzhi2012[AT]gmail[DOT].com>, BGI-Shenzhen

ENV DEBIAN_FRONTEND=noninteractive


# Install required packages
RUN apt-get update
RUN apt-get install -y  wget bzip2


# install anaconda
RUN mkdir /app

# wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \

RUN if [ ! -d /app/anaconda ]; then \
        wget -c https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/Miniconda-latest-Linux-x86_64.sh \
        -O /app/anaconda.sh && \
        bash /app/anaconda.sh -b -p /app/anaconda && \
        rm -rf /app/anaconda.sh ; fi

RUN apt-get clean && \
    apt-get autoclean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# set anaconda path
ENV PATH="/app/anaconda/bin:$PATH"


# install dependency for MitoZ
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

RUN conda install -y libgd=2.2.4 python=3.6.0 biopython=1.69 ete3=3.0.0b35  perl-list-moreutils perl-params-validate perl-clone circos=0.69 perl-bioperl blast=2.2.31  hmmer=3.1b2  bwa=0.7.12 samtools=1.3.1 infernal=1.1.1 tbl2asn openjdk ; \
    conda clean -y -a


# download MitoZ and install
RUN mkdir /mitoz_tmp && cd /mitoz_tmp && wget -c https://raw.githubusercontent.com/linzhi2013/MitoZ/master/version_2.4-alpha/release_MitoZ_v2.4-alpha.tar.bz2 &&  tar -jxvf release_MitoZ_v2.4-alpha.tar.bz2  && mv release_MitoZ_v2.4-alpha /app && rm -rf /mitoz_tmp


# install NCBI taxonomy database
RUN cd / ;  python3 -c 'from ete3 import NCBITaxa; ncbi = NCBITaxa()' ; rm -rf taxdump.tar.gz

ENV LC_ALL=C
ENV PATH=/app/anaconda/bin:$PATH

VOLUME /project
WORKDIR /project

# ENTRYPOINT ["/app/anaconda/bin/python3 /app/release_MitoZ_v2.4-alpha/MitoZ.py"]
CMD ["/bin/bash"]
