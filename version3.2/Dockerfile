FROM ubuntu:16.04

MAINTAINER Guanliang Meng <linzhi2012[AT]gmail[DOT].com>, ZFMK

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

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda install mamba -n base -c conda-forge 

# mamba install -c bioconda mitoz
# or local installation:

mamba install -f https://objects.githubusercontent.com/github-production-release-asset-2e65be/158059318/d9dc5a03-f16d-452f-9434-ec8f62df9868?X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=AKIAIWNJYAX4CSVEH53A%2F20220513%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Date=20220513T110834Z&X-Amz-Expires=300&X-Amz-Signature=676864c87db47afe694e2f4ce258812a49e3c91853058b2dafc29ec88f731dfe&X-Amz-SignedHeaders=host&actor_id=6203542&key_id=0&repo_id=158059318&response-content-disposition=attachment%3B%20filename%3DmitozEnv.yaml&response-content-type=application%2Foctet-stream

pip install https://github.com/linzhi2013/MitoZ/releases/download/3.2/mitoz-3.2.tar.gz    


# install NCBI taxonomy database
RUN cd / ;  python3 -c 'from ete3 import NCBITaxa; ncbi = NCBITaxa()' ; rm -rf taxdump.tar.gz

ENV LC_ALL=C
ENV PATH=/app/anaconda/bin:$PATH

VOLUME /project
WORKDIR /project

CMD ["/bin/bash"]
