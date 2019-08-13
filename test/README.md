This folder contains a small dataset for testing if your MitoZ is fine.

The dataset has been tested with the Docker image `guanliangmeng/mitoz:2.3`.

# How to run

If you are in the Docker group, run:

    fq1=test.1.fq.gz
    fq2=test.2.fq.gz
    outprefix=test
    docker run -v $PWD:/project  --rm guanliangmeng/mitoz:2.3 \
    python3 /app/release_MitoZ_v2.3/MitoZ.py all2 \
    --genetic_code 5 \
    --clade Arthropoda \
    --insert_size 250 \
    --thread_number 4 \
    --fastq1 $fq1 \
    --fastq2 $fq2 \
    --outprefix $outprefix \
    --fastq_read_length 125 \
    1>m.log 2>m.err


Otherwise, you may add the `sudo` command before `docker`, i.e.:


    fq1=test.1.fq.gz
    fq2=test.2.fq.gz
    outprefix=test
    sudo docker run -v $PWD:/project  --rm guanliangmeng/mitoz:2.3 \
    python3 /app/release_MitoZ_v2.3/MitoZ.py all2 \
    --genetic_code 5 \
    --clade Arthropoda \
    --insert_size 250 \
    --thread_number 4 \
    --fastq1 $fq1 \
    --fastq2 $fq2 \
    --outprefix $outprefix \
    --fastq_read_length 125 \
    1>m.log 2>m.err


# Resource usage

* 4 CPUs

* About 2G memory

* 15 minutes (wall time)


An example result is in `test.result.tar.gz`.

