#!/usr/bin/python3
"""
pick_seq_from_scoing_result_base_arr_v3.py

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

usage = """

python3 %s <nhmmer.besthit.sim.reformat.sorted> <outfile> <outfile2> <outfile3>

1. will only look for the lacking CDS for the first sequence!
2. sequences which have >= 5 full-length CDS will also be output.

""" % sys.argv[0]

if len(sys.argv) != 5:
	print(usage)
	sys.exit(0)

in_f, out_f, out_f2, out_f3 = sys.argv[1:5]

gene_order = ['COX1', 'COX2', 'COX3', 'ND1', 'ND2', 'ND3','ND4', 'ND4L', 'ND5', 'ND6', 'CYTB', 'ATP6', 'ATP8']
gene_dict = {'COX1':0, 'COX2':1, 'COX3':2, 'ND1':3, 'ND2':4, 'ND3':5,
       		'ND4':6, 'ND4L':7, 'ND5':8, 'ND6':9, 'CYTB':10, 'ATP6':11, 'ATP8':12}


def gene2arr(seq):
	edge_genes = []
	internal_genes = []
	for line in seq:
		i = line.split("\t")
		gene = i[0]
		gene_length_enough = int(i[-1])
		if gene_length_enough == 4:
			internal_genes.append(gene)
		elif gene_length_enough == 2:
			edge_genes.append(gene)
		else:
			pass
	gene_13 = ['.' for i in range(13)]
	for i in edge_genes:
		j = gene_dict[i]
		gene_13[j] = 2
	for i in internal_genes:
		j = gene_dict[i]
		gene_13[j] = 4

	#print(gene_13)
	return gene_13


######## readin scoring result ########
rec_dict = {}
gene_two_dim = []

fh_in2 = open(in_f, 'r')
firstline = 1
idline = ''
seq = []
id_order = []
for i in fh_in2:
	i = i.strip()
	#print(i)
	if i.startswith(">"):
		if firstline:
			firstline = 0
		else:
			rec_dict[idline] = seq
			gene_two_dim.append(gene2arr(seq))
			id_order.append(idline)
		idline = i
		seq = []
	else:
		seq.append(i)

rec_dict[idline] = seq
gene_two_dim.append(gene2arr(seq))
id_order.append(idline)

fh_in2.close()

#############

seq_hmm_dict = {}

for seqid in id_order:
	for line in rec_dict[seqid]:
		i = line.split("\t")
		gene = i[0]
		hmm_start, hmm_end = i[3:5]
		hmm_start = int(hmm_start)
		hmm_end = int(hmm_end)
		if seqid not in seq_hmm_dict:
			seq_hmm_dict.setdefault(seqid, {})
		if gene not in seq_hmm_dict[seqid]:
			seq_hmm_dict[seqid].setdefault(gene, [])
		seq_hmm_dict[seqid][gene].extend([hmm_start, hmm_end])



##################### reformat the nhmmer result ################

fh_out = open(out_f, 'w')

for gene in gene_order:
	print(gene, end="\t", file=fh_out)
print(file=fh_out)

for i in range(0, len(gene_two_dim)):
	for j in gene_two_dim[i]:
		print(j, end="\t", file=fh_out)
	print(id_order[i], file=fh_out)
fh_out.close()
########### determine which sequences should be output ##########


#gene_order = ['COX1', 'COX2', 'COX3', 'ND1', 'ND2', 'ND3','ND4', 'ND4L', 'ND5', 'ND6', 'CYTB', 'ATP6', 'ATP8']
gene_13_hmmPostion_dict = {}
for gene in gene_order:
	if gene not in gene_13_hmmPostion_dict:
		gene_13_hmmPostion_dict.setdefault(gene, [0, 0])

# gene_dict = {'COX1':0, 'COX2':1, 'COX3':2, 'ND1':3, 'ND2':4, 'ND3':5,
#       		'ND4':6, 'ND4L':7, 'ND5':8, 'ND6':9, 'CYTB':10, 'ATP6':11, 'ATP8':12}

id_output = []
gene_13 = [0 for i in range(13)]
firstline = 1
for i in range(0, len(gene_two_dim)):
	seqid = id_order[i]
	keep = 1
	if firstline:
		firstline = 0
		id_output.append(id_order[i])
		for j in range(13):
			gene = gene_order[j]
			if gene_two_dim[i][j] != '.':
				gene_13[j] = gene_two_dim[i][j]
				gene_13_hmmPostion_dict[gene] = seq_hmm_dict[seqid][gene]
	else:
		tmp = 0
		count1 = 0
		count2 = 0
		for j in range(13):
			gene = gene_order[j]
			if (gene_13[j] == 0) and (gene_two_dim[i][j] in [2, 4]):
				count1 += 1
			elif (gene_13[j] == 4 and gene_two_dim[i][j] in [2, 4]) or (gene_13[j] == 2 and gene_two_dim[i][j] == 4):
				count2 += 1
			elif (gene_13[j] == gene_two_dim[i][j] == 2):
				## according to hmm_start and hmm_end
				condiction1 = (gene_13_hmmPostion_dict[gene][0] + 30) <= seq_hmm_dict[seqid][gene][0] <= (gene_13_hmmPostion_dict[gene][1] - 30)
				condiction2 = (gene_13_hmmPostion_dict[gene][0] + 30) <= seq_hmm_dict[seqid][gene][1] <= (gene_13_hmmPostion_dict[gene][1] - 30)
				if condiction1 or condiction2:
					count2 += 1

		## brach 2, output if >= 5, but not update gene_13
		if count2 > 0:
			if gene_two_dim[i].count(4) >= 5:
				#id_output.append(id_order[i]) # original
				id_output.append(id_order[i]+' why=OutputBecauseOfNoLessThanFivePCGs')
			continue
		## brach 1, output, and update gene_13
		elif count1 > 0:
			id_output.append(id_order[i])
			for j in range(13):
				gene = gene_order[j]
				if (gene_13[j] == 0) and (gene_two_dim[i][j] in [2, 4]):
					gene_13[j] = gene_two_dim[i][j]
					gene_13_hmmPostion_dict[gene] = seq_hmm_dict[seqid][gene]
				elif gene_13[j] == gene_two_dim[i][j] == 2:
					gene_13[j] += 2
					gene_13_hmmPostion_dict[gene] = seq_hmm_dict[seqid][gene]
			continue

fh_out2 = open(out_f2, 'w')
for i in id_output:
	print(i, file=fh_out2)
	if ' why=OutputBecauseOfNoLessThanFivePCGs' in i:
		i = i.replace(' why=OutputBecauseOfNoLessThanFivePCGs', '')
	for j in rec_dict[i]:
		print(j, file=fh_out2)
fh_out2.close()


## output genes found statistic
fh_out3 = open(out_f3, 'w')

print("                ", end="\t", file=fh_out3)
for i in range(13):
	#print("full length CDS:",i+1, end="\t")
	print(i+1, sep="", end="\t", file=fh_out3)
print(file=fh_out3)

full_len_count = gene_13.count(4)
print(full_len_count, "full length CDS:", end="\t", file=fh_out3)
for j in range(13):
	if gene_13[j] == 4:
		print(gene_order[j], end="\t", file=fh_out3)
print(file=fh_out3)

half_len_count = gene_13.count(2)
print(half_len_count, "half length CDS:", end="\t", file=fh_out3)
for j in range(13):
	if gene_13[j] == 2:
		print(gene_order[j], end="\t", file=fh_out3)
print(file=fh_out3)

not_fount_count = gene_13.count(0)
print(not_fount_count, "CDS not found:", end="\t", file=fh_out3)
for j in range(13):
	if gene_13[j] == 0:
		print(gene_order[j], end="\t", file=fh_out3)
print(file=fh_out3)










