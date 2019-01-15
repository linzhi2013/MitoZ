#!/usr/bin/python3
"""
reformat_nhmmer_besthit-sim.py

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
import re

usage = """

python3 %s <test.hmmout2tbl.besthit.sim> <outfile>

output format:
gene, seq_start, seq_end, hmm_start, hmm_end, strand

""" % sys.argv[0]

if len(sys.argv) != 3:
    print(usage)
    sys.exit(0)

infile, outfile = sys.argv[1:3]

seq_dict = {}

abundance_dict = {}

fh_in = open(infile, 'r')
for i in fh_in:
    i = i.strip()
    i = i.split("\t")
    seqid, gene, seq_desc = [i[k] for k in (0, 1,  -1)]
    seq_len = "length=" + i[6]
    seqid = seqid + " " + seq_desc + " " + seq_len

#   if seqid.startswith("C"):
#       abundance = float(i[-1])
#   else:
#       abundance = i[-1].split()[-2]
#       abundance = float(abundance)
    abundance = re.search(r'(\d+[\.]*\d+)', i[-1]).group(1)
    abundance_dict[seqid] = abundance
    # gene, seq_start, seq_end, hmm_start, hmm_end, strand
    j = [i[k] for k in (1, 4, 5, 2, 3, 7)]
    j[1] = int(j[1])
    j[2] = int(j[2])
    if j[1] > j[2]:
        j[1], j[2] = j[2], j[1]

    seq_start = j[1]
    j = [str(k) for k in j]
    j = "\t".join(j)

    if seqid not in seq_dict:
        seq_dict.setdefault(seqid, {})
    if gene not in seq_dict[seqid]:
        seq_dict[seqid].setdefault(gene, [])

    seq_dict[seqid][gene].extend([seq_start, j])

fh_in.close()

seqid_sorted = []
for seqid, abundance in sorted(abundance_dict.items(), key=lambda x: x[1], reverse=True):
    seqid_sorted.append(seqid)

fh_out = open(outfile, 'w')
for seqid in seqid_sorted:
    print(">%s" %seqid, file=fh_out)
    for gene, (seq_start, j) in sorted(seq_dict[seqid].items(), key=lambda x: x[1][0]):
        print(j, file=fh_out)
fh_out.close()








