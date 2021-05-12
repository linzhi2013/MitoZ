#!/usr/bin/env python3
import sys
import os
import re
import collections
from Bio import SeqIO


def get_fasta_stat(fa_file=None, most_related_sp_file=None):
    seqid_len_topology_relatedSP = collections.defaultdict(list)
    for rec in SeqIO.parse(fa_file, 'fasta'):
        seqid = str(rec.id)
        seqlen = str(len(rec))
        toplogy = 'no'
        if 'topology=circular' in rec.description:
            toplogy = 'yes'

        seqid_len_topology_relatedSP[seqid].extend([seqlen, toplogy])

    with open(most_related_sp_file, 'r') as fh:
        for i in fh:
            i = i.strip()
            seqid, sp = i.split('\t')
            seqid_len_topology_relatedSP[seqid].append(sp)

    # print(seqid_len_topology_relatedSP)

    return seqid_len_topology_relatedSP


def main():
    usage = '''
python3 {0} <in.fa> <most_related_sp_file>
'''.format(sys.argv[0])

    if len(sys.argv) != 3:
        sys.exit(usage)

    fa_file, most_related_sp_file = sys.argv[1:3]

    seqid_len_topology_relatedSP = get_fasta_stat(fa_file=fa_file, most_related_sp_file=most_related_sp_file)


    col_fomrat = '{:15s}{:15s}{:15s}{:20s}'
    tpo_title_line = col_fomrat.format('#Seq_id', 'Length(bp)', 'Circularity', 'Closely_related_species')
    print(tpo_title_line)
    for seqid in seqid_len_topology_relatedSP:
        seqlen, toplogy, sp = seqid_len_topology_relatedSP[seqid]
        out = col_fomrat.format(seqid, seqlen, toplogy, sp)
        print(out)

    print("\nIn the next step, you may want to run the 'annotate' module of MitoZ.\nIn case there are missing genes, you might try to find them from the\n'*.high_abundance*' and '*.low_abundance*' files!")



if __name__ == '__main__':
    main()
