fq1=test.1.fq.gz
fq2=test.2.fq.gz
outprefix=test
time docker run -v $PWD:/project  --rm guanliangmeng/mitoz:2.3 \
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
