#!/usr/bin/python3
"""
sort_gb_by_len.py

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
import argparse

desc = """
sort sequences in Genbank/FASTA files by sequence length.
"""

parser = argparse.ArgumentParser(description=desc)

parser.add_argument("-i", metavar="<STR>", required=True,
					help="in Genbank/FASTA file")
parser.add_argument("-o", metavar="<STR>", required=True,
					help="out Genbank/FASTA file")
parser.add_argument("-l", default=False, action="store_true",
					help="sort by decrease length (default: sort by increase)")
parser.add_argument("-f", default=False, action="store_true",
					help="input is FASTA file (default: Genbank file)")

if len(sys.argv) == 1:
	parser.print_help()
	parser.exit()
else:
	args = parser.parse_args()


rec_dict = {}
seqlen_dict = {}

seq_format = 'genbank'
if args.f:
	seq_format = 'fasta'

## to avoid same seqids in records
j = 1
for rec in SeqIO.parse(args.i, seq_format):
	seqid = "seq" + str(j)
	j += 1
	seqlen_dict[seqid] = len(rec.seq)
	rec_dict[seqid] = rec

if args.l:
	ids_sorted = sorted(seqlen_dict.items(), key=lambda x:x[1], reverse=True)
else:
	ids_sorted = sorted(seqlen_dict.items(), key=lambda x:x[1])

print(ids_sorted)

fh_out = open(args.o, 'w')
for (seqid, seqlen) in ids_sorted:
	SeqIO.write(rec_dict[seqid], fh_out, seq_format)

fh_out.close()
