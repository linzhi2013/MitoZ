#!/usr/bin/python3
"""
filter_mitobim-bases_by_depth.py

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

filter out bases from MitoBim which have <= cutoff-1 sequencing depth.

output result sequences composition: 5' part from mitobim + SOAPTrans + 3' part from mitobim


python3 %s  <SOAPTrans-MitoBim merged sequence file (fasta)>  <sequencing depth file (one-line per sequence)>  <outfile>  [sequencing depth cutoff]

sequencing depth file format    seqid X1 X2 ... Xn
sequencing depth cutoff    default: 2, means all output bases have a sequencing depth >= 2X.

""" % sys.argv[0]

if len(sys.argv) not in [4, 5]:
	print(usage)
	sys.exit(0)

depth_cutoff = 2

merged_file, depth_file, outfile = sys.argv[1:4]
if len(sys.argv) == 5:
	depth_cutoff = int(sys.argv[4])

seq_dict = {}
soaptrans_pos = {}

depth_dict = {}


## 考虑merged序列组成：
# 1. mitobim + SOAPtrans + mitobim
# 2. SOAPtrans + mitobim
# 3. mitobim + SOAPtrans

for rec in SeqIO.parse(merged_file, 'fasta'):
	seq = str(rec.seq)
	seq_dict[rec.id] = seq
	soaptrans_pos.setdefault(rec.id, [-1, -1])

	## 1. mitobim + SOAPtrans + mitobim
	tmp = 0
	if seq[0].islower() and seq[-1].islower():
		for i in range(0, len(seq)):
			if (tmp == 0) and seq[i].isupper():
				tmp = 1
				soaptrans_pos[rec.id][0] = i
			if tmp == 1 and seq[i+1].islower():
				soaptrans_pos[rec.id][1] = i
				break

	## 2. SOAPtrans + mitobim
	if seq[0].isupper() and seq[-1].islower():
		for i in range(0, len(seq)):
			if seq[i+1].islower():
				soaptrans_pos[rec.id][1] = i
				break

	## 3. mitobim + SOAPtrans
	if seq[0].islower() and seq[-1].isupper():
		for i in range(0, len(seq)):
			if seq[i].isupper():
				soaptrans_pos[rec.id][0] = i
				break


fh_in = open(depth_file, 'r')
for i in fh_in:
	i = i.strip().split()
	seqid = i[0]
	j = [int(k) for k in i[1:]]
	if seqid not in depth_dict:
		depth_dict.setdefault(seqid, [])
	depth_dict[seqid].extend(j)
fh_in.close()


fh_out = open(outfile, 'w')
for seqid in soaptrans_pos.keys():
	mitobim_5_end_min_pos = -1
	mitobim_3_end_min_pos = -1

	## 1. mitobim + SOAPtrans + mitobim
	if 0 < soaptrans_pos[seqid][0] < soaptrans_pos[seqid][1]:
		trans_start = soaptrans_pos[seqid][0]
		trans_end = soaptrans_pos[seqid][1]
		mitobim_5_end_min = min(depth_dict[seqid][0:trans_start])
		mitobim_3_end_min = min(depth_dict[seqid][trans_end:])

		#print(mitobim_5_end_min, mitobim_3_end_min)

		## if the minimum depth less than depth_cutoff
		if mitobim_5_end_min < depth_cutoff:
			for i in range(0,trans_start):
				if depth_dict[seqid][i] < depth_cutoff:  # Note, here should be depth_cutoff, but not mitobim_5_end_min
					mitobim_5_end_min_pos = i  # the last position
		if mitobim_3_end_min < depth_cutoff:
			for i in range(trans_end, len(depth_dict[seqid])):
				if depth_dict[seqid][i] < depth_cutoff:
					mitobim_3_end_min_pos = i
					break  # the first position

	## 2. SOAPtrans + mitobim
	if soaptrans_pos[seqid][0] < 0 < soaptrans_pos[seqid][1]:
		trans_end = soaptrans_pos[seqid][1]
		mitobim_3_end_min = min(depth_dict[seqid][trans_end:])
		if mitobim_3_end_min < depth_cutoff:
			for i in range(trans_end, len(depth_dict[seqid])):
				if depth_dict[seqid][i] < depth_cutoff:
					mitobim_3_end_min_pos = i
					break  # the first position

	## 3. mitobim + SOAPtrans
	if soaptrans_pos[seqid][1] < 0 < soaptrans_pos[seqid][0]:
		trans_start = soaptrans_pos[seqid][0]
		mitobim_5_end_min = min(depth_dict[seqid][0:trans_start])
		if mitobim_5_end_min < depth_cutoff:
			for i in range(0, trans_start):
				if depth_dict[seqid][i] < depth_cutoff:
					mitobim_5_end_min_pos = i  # the last position

	## output
	if (mitobim_5_end_min_pos == -1) and (mitobim_3_end_min_pos == -1):
		print(">%s" % seqid, file=fh_out)
		print(seq_dict[seqid], file=fh_out)

	if (mitobim_5_end_min_pos == -1) and (mitobim_3_end_min_pos > 0):
		print(">%s" % seqid, file=fh_out)
		print(seq_dict[seqid][0:mitobim_3_end_min_pos], file=fh_out)

	if (mitobim_5_end_min_pos > 0) and (mitobim_3_end_min_pos == -1):
		print(">%s" % seqid, file=fh_out)
		print(seq_dict[seqid][(mitobim_5_end_min_pos+1):], file=fh_out)

	if 0 < mitobim_5_end_min_pos < mitobim_3_end_min_pos:
		print(">%s" % seqid, file=fh_out)
		print(seq_dict[seqid][(mitobim_5_end_min_pos+1):mitobim_3_end_min_pos], file=fh_out)

fh_out.close()



