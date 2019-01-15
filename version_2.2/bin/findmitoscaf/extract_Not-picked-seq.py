#!/usr/bin/env python3
import sys
import os
import re
from Bio import SeqIO


def get_not_picked_score(score_file=None, score_picked=None):
    score_picked_seqids = []
    with open(score_picked, 'r') as fh:
        for i in fh:
            i = i.strip()
            if not i.startswith('>'):
                continue
            seqid = i.split()[0].replace('>', '')
            score_picked_seqids.append(seqid)

    score_not_picked = {}
    with open(score_file, 'r') as fh:
        first = True
        seq_lines = []
        for i in fh:
            i = i.strip()
            if i.startswith('>'):
                if first:
                    first = False
                else:
                    if seqid not in score_picked_seqids:
                        score_not_picked[seqid] = seq_lines

                seqid = i.split()[0].replace('>', '')
                seq_lines = []
                seq_lines.append(i)
            else:
                seq_lines.append(i)

        # the last
        if seqid not in score_picked_seqids:
            score_not_picked[seqid] = seq_lines

    return score_not_picked


def extract_not_picked_seq(in_fasta=None, score_not_picked=None, out_not_picked_fa=None):

    fh_out_fa = open(out_not_picked_fa, 'w')
    for rec in SeqIO.parse(in_fasta, 'fasta'):
        if str(rec.id) in score_not_picked:
            SeqIO.write(rec, fh_out_fa, 'fasta')

    fh_out_fa.close()


def main():
    usage = '''
python3 {0} <*.high_abundance_10.0X.reformat.sorted> <*.high_abundance_10.0X.reformat.sorted.picked> <DRR095708.hmmtblout.besthit.sim.filtered.fa> <out.Not_picked> <out.Not_picked.fa>
'''.format(sys.argv[0])

    if len(sys.argv) != 6:
        sys.exit(usage)

    score_file, score_picked, in_fasta, out_not_picked, out_not_picked_fa = sys.argv[1:6]

    score_not_picked = get_not_picked_score(score_file=score_file, score_picked=score_picked)

    extract_not_picked_seq(in_fasta=in_fasta, score_not_picked=score_not_picked, out_not_picked_fa=out_not_picked_fa)

    fhout = open(out_not_picked, 'w')
    for seqid in score_not_picked:
        lines = '\n'.join(score_not_picked[seqid])
        print(lines, file=fhout)
    fhout.close()

if __name__ == '__main__':
    main()

