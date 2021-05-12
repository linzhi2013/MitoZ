#!/usr/bin/python3
"""
prepare_fsa_v2.py

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

    Add some necessary information into defline of fasta file.
    output sequences are sorted from minimum length to maximum length.

Usage

    python3 {0} <candidate.mito.target_taxon.fa> <organism> <circular|linear> <table> <outfile>

""".format(sys.argv[0])

if len(sys.argv) != 6:
    print(usage)
    sys.exit(0)

seq_f, organism, topology, mgcode, out_f = sys.argv[1:6]

seq_dict = {}
seqlen_dict = {}

for rec in SeqIO.parse(seq_f, 'fasta'):
	seqlen_dict[rec.id] = len(rec.seq)
	seq_dict[rec.id] = rec.seq

seqids_sorted = sorted(seqlen_dict.items(), key=lambda x:x[1])

print(seqids_sorted)

fh_out = open(out_f, 'w')
for (seqid, seqlen) in seqids_sorted:
	print(">"+seqid, "[organism=%s]" % organism,
		 "[topology=%s]" % topology, "[mgcode=%s]" % mgcode,
		 "[location=mitochondrion]\n"+seq_dict[seqid], file=fh_out)

#fh_in = open(seq_f, 'r')
#fh_out = open(out_f, 'w')
#for line in fh_in:
#    line = line.rstrip()
#    if line.startswith(">"):
#        print(line, "[organism=%s]" % organism, "[topology=%s]" % topology, "[mgcode=%s]" % mgcode, "[location=mitochondrion]", file=fh_out)
 #   else:
 #       print(line, file=fh_out)

#fh_in.close()
fh_out.close()
