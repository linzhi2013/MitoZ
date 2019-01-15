#!/usr/bin/python3
"""
circle_check.py

Copyright (c) 2017-2018 Guanliang Meng <mengguanliang@foxmail.com>.

This file is part of MitoZ.

MitoZ is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MitoZ is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MitoZ.  If not, see <http://www.gnu.org/licenses/>.

"""

import sys
from Bio import SeqIO

usage = """
Description

    Checking whether the sequences are circular when the sequences have
    length >= 12Kbp

Usage

    python3 {0}  <in.fasta>  <outPrefix> <mismatch_allowed>

output files:

1. <outPrefix>.mitogenome.fa
All the sequences from <in.fasta>.

The sequence id line will be like:
>C1 topology=circular
>C2 topology=linear

For the circular mt sequence, the overlapping region (the second `ATGCNN`
below) has been removed (below is an example)

ATGCNNNNN[ATGCNN]

Assuming `ATGCNNNNN` is a circular mt sequence, `ATGCNN` are the overlapping
regions.

2. <outPrefix>.start2end_for-circular-mt-only
This file contains the circular sequences only, and the first 300 bp of each
has been moved to the end of the sequence, just for better reads mapping. You
can check the sequencing depth around the 'joining site' (-300 bp) using the
`annotate` module of MitoZ, to confirm if the sequence is really circular.

3. <outPrefix>.overlap_information
The overlapping sequence detected for the circular sequences.

""".format(sys.argv[0])

if len(sys.argv) != 4:
    print(usage)
    sys.exit(0)

in_f, outPrefix, mismatch_allowed = sys.argv[1:4]
mismatch_allowed = int(mismatch_allowed)

def seq_distance(seq1=None, seq2=None):
    if len(seq1) != len(seq2):
        print(seq1)
        print(seq2)
        sys.exit("Error: seq1 and seq2 have different lengths!")

    seq1 = seq1.upper()
    seq2 = seq2.upper()

    mismatch = 0
    for i in range(0, len(seq1)):
        if seq1[i] != seq2[i]:
            mismatch += 1
    return (mismatch)

fh_out = open(outPrefix+'.mitogenome.fa', 'w')
fh_out2 = open(outPrefix+'.overlap_information', 'w')
fh_out3 = open(outPrefix+".start2end_for-circular-mt-only", 'w')

for rec in SeqIO.parse(in_f, 'fasta'):
    seq = str(rec.seq)
    seqlen = len(seq)
    if seqlen < 12000:
        outline = ">" + rec.id + ' topology=linear'+ "\n" + seq
        print(outline, file=fh_out)
        continue

    potential_seq = seq[12000:]
    overlap_len = 5
    match_start_pos = 0
    tmp = 0
    while(1):
        subseq_5 = seq[0:overlap_len]
        match_count = potential_seq.count(subseq_5)
        if match_count > 1:
            overlap_len += 1
        elif match_count == 1:
            match_start_pos = 12000 + potential_seq.find(subseq_5)
            #print("first: ", subseq_5)
            break
        else:
            print("not found!")
            tmp = 1
            break

    if tmp == 1:
        outline = ">" + rec.id + ' topology=linear'+ "\n" + seq
        print(outline, file=fh_out)
        continue

    start_pos = match_start_pos
    overlap_seq = ""
    for i in range(0, int(seqlen/2)):
        subseq_5 = seq[0:overlap_len]
        stop_pos = match_start_pos + overlap_len
        if stop_pos <= len(seq):
            subseq_3 = seq[start_pos:stop_pos]
            if seq_distance(subseq_5, subseq_3) < mismatch_allowed:
                overlap_seq = subseq_5
                overlap_len += 1
            else:
                print("A: not found!")
                outline = ">" + rec.id + ' topology=linear'+ "\n" + seq
                print(outline, file=fh_out)
                break ## no circular found, go to next sequence
        else:
            outline2 = ">" + rec.id + " overlap between 5' and 3' are " + \
                    str(overlap_len) +"bp\n" + overlap_seq
            print(outline2, file=fh_out2)

            outline = ">" + rec.id + ' topology=circular'+ "\n" + seq[0:start_pos]
            print(outline, file=fh_out)

            outline3 = ">" + rec.id + ' topology=circular' + "\n" + seq[300:start_pos]+seq[0:300]
            print(outline3, file=fh_out3)

            break

fh_out.close()

