#!/usr/bin/python3
"""
pick_target_taxon_seq.py

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

    To extract sequnces from 'candidate.mito.fa' based on 'candidate.CDS.position.blastn_nt.besthit.sim' and 'ct5.gi.taxon.Arthropoda'.

Usage

    python3 {0}  <ct5.gi.taxon.Arthropoda> <candidate.mito.nt.besthit.sim>  <candidate.mito.fa>  <outfile>

""".format(sys.argv[0])

if len(sys.argv) != 5:
    print(usage)
    sys.exit(0)

taxon_f, blast_f, seq_f, out_f = sys.argv[1:5]


need_gis = []
need_scafs = {}

fh_in1 = open(taxon_f, 'r')
for line in fh_in1:
    line = line.rstrip()
    gi = line.split("\t")[0]
    need_gis.append(gi)

fh_in1.close()

fh_in2 = open(blast_f, 'r')
for line in fh_in2:
    line = line.rstrip()
    query, gi = line.split("\t")[0:2]
    if gi in need_gis:
        need_scafs[query] = 1

fh_in2.close()

count = len(need_scafs)

fh_out = open(out_f, 'w')
for rec in SeqIO.parse(seq_f, 'fasta'):
    if count < 1:
        break

    if rec.id in need_scafs:
        SeqIO.write(rec, fh_out, 'fasta')
        need_scafs[rec.id] = 0
        count -= 1

fh_out.close()

for seq_id in need_scafs.keys():
    if need_scafs[seq_id] == 1:
        print("can not find %s of %s in %s" % (seq_id, blast_f, seq_f))
