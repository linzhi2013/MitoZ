#!/usr/bin/python3
import collections
"""
extract_nhmmer_seq.py

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

import argparse
from Bio import SeqIO

description = "extract fasta sequnces of simplified nhmmer file"

parser = argparse.ArgumentParser(description=description)

parser.add_argument("-f", required=True, help="simplified nhmmer result " + \
					" file (.tbl format)")
parser.add_argument("-d", required=True, help="input fasta file")
parser.add_argument("-o", required=True, help="outfile")

args = parser.parse_args()

#seqid_dict = {}
seqid_dict = collections.defaultdict(list)

fh_in = open(args.f, 'r')
for i in fh_in:
	i = i.rstrip()
	if not i.startswith(">"):
		continue
	line = i.split()
	seqid = line[0].replace(">", "")
	others = ""
	if len(line) > 1:
		for j in line[1:]:
			# length=6397 score=422.695 why=
			if 'length=' in j or 'score=' in j or 'why=' in j:
				others += ' ' + j
	#print(others)
	seqid_dict[seqid] = [1, others]

fh_in.close()

# print(seqid_dict)

count = len(seqid_dict)
fh_out = open(args.o, 'w')
for rec in SeqIO.parse(args.d, 'fasta'):
	if count <= 0:
		break
	if rec.id in seqid_dict:
		seqid_dict[rec.id][0] = 0
		count -= 1

		others = seqid_dict[rec.id][1]
		rec.description += others
		if 'why=' in others:
			rec.id += '_FivePCGs'
		SeqIO.write(rec, fh_out, "fasta")
		

fh_out.close()


for seqid in seqid_dict.keys():
	if seqid_dict[seqid][0] == 1:
		print("can not find seqid %s in file %s" % (seqid, args.d))
