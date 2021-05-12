#!/usr/bin/python3
"""
scoring_nhmmer_besthit_sim_reformat.py

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

python3 %s  <CDS_length_list> <nhmmer.besthit.sim.reformat> <outfile>

caculate the average CDS annotated length rate, then output the sorted result.

""" % sys.argv[0]

if len(sys.argv) != 4:
	print(usage)
	sys.exit(0)

CDS_length_list, nhmmer_file, outfile = sys.argv[1:4]

defined_CDS_length_dict = {}

fh_in = open(CDS_length_list, 'r')
for i in fh_in:
	i = i.strip()
	gene, gene_len = i.split()
	defined_CDS_length_dict[gene] = int(gene_len)
fh_in.close()

#################

def cal_score(idline, seq):
	abundance = 0.0
	if idline.startswith(">C"):
		abundance = idline.split()[1]
	else:
		abundance = idline.split()[2]
	abundance = float(abundance)

	seq_len = idline.split("length=")[1]
	seq_len = int(seq_len)


	###### total score #######
	total_score = 0.0
	total_annotated_rate = 0.0

	for line in seq:
		gene, seq_start, seq_end, hmm_start, hmm_end, strand = line.split("\t")
		#gene = re.split(r'[.-]', gene)[-1]
		seq_start, seq_end = [int(i) for i in (seq_start, seq_end)]
		annotated_len = abs(seq_start - seq_end + 1)
		annotated_rate = annotated_len / defined_CDS_length_dict[gene]
		total_annotated_rate += annotated_rate

	total_score = total_annotated_rate * abundance
	#total_score = round(total_score, 3)

	####### total_length_enough ## start ############

	total_length_enough = [0 for k in range(len(seq))]

	def start_OR_end_gene_length_enough(idline, line):
		gene, seq_start, seq_end, hmm_start, hmm_end, strand = line.split("\t")
		#gene = re.split(r'[.-]', gene)[-1]
		seq_start, seq_end, hmm_start, hmm_end = [int(i) for i in (seq_start, seq_end, hmm_start, hmm_end)]
		annotated_len = abs(seq_start - seq_end + 1)

		## 4: full_length, 2: not_full_length, 0: unannotated
		one_length_enough = 2

		if annotated_len >= defined_CDS_length_dict[gene]:
			one_length_enough = 4
		else:
			if strand == '+':
				seq_start = seq_start - (defined_CDS_length_dict[gene] - annotated_len)
				if seq_start > 0:
					one_length_enough = 4
			elif strand == '-':
				seq_end = seq_end + (defined_CDS_length_dict[gene] - annotated_len)
				if seq_end <= seq_len:
					one_length_enough = 4

		return one_length_enough

	if len(seq) == 1:
		# one gene
		total_length_enough[0] = start_OR_end_gene_length_enough(idline, seq[0])
	elif len(seq) >= 2:
		## start gene
		total_length_enough[0] = start_OR_end_gene_length_enough(idline, seq[0])
		## end gene
		total_length_enough[-1] = start_OR_end_gene_length_enough(idline, seq[-1])
		## internal genes
		for k in range(1, len(seq)):
			total_length_enough[k] = 4

	####### total_length_enough ## end ############

	return total_score, abundance, total_length_enough


#################
total_score_dict = {}
rec_dict = {}

fh_in2 = open(nhmmer_file, 'r')
firstline = 1
for i in fh_in2:
	i = i.strip()
	if i.startswith(">"):
		if firstline:
			firstline = 0
		else:
			total_score_dict[idline] = cal_score(idline, seq)
			rec_dict[idline] = seq
		idline = i
		seq = []
	else:
		seq.append(i)

total_score_dict[idline] = cal_score(idline, seq)
rec_dict[idline] = seq

fh_in2.close()

## idline, [total_score, abundance, total_length_enough]
sorted_by_abundance = sorted(total_score_dict.items(), key=lambda x: x[1][1], reverse=True)

fh_out = open(outfile, 'w')
for a, b in sorted(sorted_by_abundance, key=lambda x: x[1][0], reverse=True):
	print(a, "score=%.3f" % b[0], file=fh_out)
	k = 0
	for line in rec_dict[a]:
		one_length_enough = total_score_dict[a][2][k]
		k += 1
		print(line, one_length_enough, sep="\t", file=fh_out)
fh_out.close()
