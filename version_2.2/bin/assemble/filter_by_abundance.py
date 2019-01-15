#!/usr/bin/env python3
import sys
import os
import re
from Bio import SeqIO


def filter_by_abundance(in_sim=None, in_fasta=None, min_abundance=10):
    prefix = os.path.basename(in_sim)
    high_abun_sim = '{0}.high_abundance_{1}X'.format(prefix, min_abundance)
    low_abun_sim = '{0}.low_abundance'.format(prefix)
    low_abun_seqids = set()
    low_abun_fasta = '{0}.low_abundance.fasta'.format(prefix)
    fhout_high = open(high_abun_sim, 'w')
    fhout_low = open(low_abun_sim, 'w')
    with open(in_sim, 'r') as fh:
        for i in fh:
            i = i.strip()
            if not i:
                continue
            line = i.split('\t')
            seqid = line[0]
            m = re.search(r'(\d+[\.]*\d+)', line[-1])
            abun = m.group(1)
            #if 'Locus' in abun or len(abun.split()) > 1:
            #    abun = abun.split()[1]
            abun = float(abun)
            if abun < min_abundance:
                print(i, file=fhout_low)
                low_abun_seqids.add(seqid)
            else:
                print(i, file=fhout_high)
    fhout_high.close()
    fhout_low.close()

    if len(low_abun_seqids) > 0:
        fhout = open(low_abun_fasta, 'w')
        for rec in SeqIO.parse(in_fasta, 'fasta'):
            if str(rec.id) in low_abun_seqids:
                SeqIO.write(rec, fhout, 'fasta')
        fhout.close()


def main():
    usage = '''
python3 {0} <min_abundance> <work71.hmmtblout.besthit.sim[.filtered]> <in.fasta>
'''.format(sys.argv[0])

    if len(sys.argv) != 4:
        sys.exit(usage)

    min_abundance, in_sim, in_fasta = sys.argv[1:4]
    filter_by_abundance(
        in_sim=in_sim,
        in_fasta=in_fasta,
        min_abundance=float(min_abundance)
    )


if __name__ == '__main__':
    main()











