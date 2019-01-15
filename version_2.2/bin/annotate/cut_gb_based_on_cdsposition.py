#!/usr/bin/python3
"""
cut_gb_based_on_cdsposition.py

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
import argparse
from Bio import SeqIO
import re

description = """
1. select the sequences have duplicate gene names in cds.position.sorted.revised, assuming they are cuicular sequences.
 2. cut these genbank records based on coordinances of the duplicate genes.
 3. non-circular genbank records will be output originally.
"""

parser = argparse.ArgumentParser(description=description)

parser.add_argument("-cps", metavar="<FILE>",
					help="cds.position.sorted.revised file")
parser.add_argument("-gb", metavar="<FILE>", help="genbank file")
parser.add_argument("-o", metavar="<FILE>", help="output file")

if len(sys.argv) == 1:
	parser.print_help()
	parser.exit()
else:
	args = parser.parse_args()

all_seq = {}
circular_seq = {}

with open(args.cps, 'r') as fh_in:
	firstline = 1
	oldseq = 0
	seqid = ""
	for i in fh_in:
		i = i.strip()
		if i.startswith(">"):
			if firstline:
				firstline = 0

			seqid = i.replace(">", "")
			all_seq.setdefault(seqid, {})
			oldseq = 0
		else:
			if oldseq:
				continue

			gene, start, end = i.split("\t")[0:3]
			start = re.sub(r">|<", "", start)
			end = re.sub(r">|<", "", end)
			start = int(start)
			end = int(end)
			gene = gene.split("-")[0]

			if gene in all_seq[seqid]:
				gene_len1 = \
					abs(all_seq[seqid][gene][0] - all_seq[seqid][gene][1])
				gene_len2 = abs(start-end)

				if gene_len1 > gene_len2:
					circular_start = all_seq[seqid][gene][0]
					circular_end = start
				else:
					circular_start = all_seq[seqid][gene][1]
					circular_end = end

				circular_seq.setdefault(seqid, [])
				circular_seq[seqid] = [circular_start, circular_end]
				#print(circular_start, circular_end)
				oldseq = 1
			else:
				all_seq[seqid].setdefault(gene, [])
				all_seq[seqid][gene].extend([start, end])

fh_out = open(args.o, 'w')
for rec in SeqIO.parse(args.gb, 'gb'):
	if (rec.id in all_seq) and (rec.id not in circular_seq):
		SeqIO.write(rec, fh_out, 'gb')
		continue
	if rec.id in circular_seq:
		circular_start, circular_end = circular_seq[rec.id]
		print(rec.id, "extract region: ", circular_start, circular_end)
		rec2 = rec[circular_start-1:circular_end]
		rec2.annotations["source"] = rec.annotations["source"]
		rec2.annotations["organism"] = rec.annotations["organism"]
		SeqIO.write(rec2, fh_out, 'gb')
fh_out.close()
