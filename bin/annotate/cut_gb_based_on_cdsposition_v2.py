#!/usr/bin/python3
"""
cut_gb_based_on_cdsposition_v2.py

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
 2. cut these genbank records based on coordinances of the duplicate
 	genes (gene distance > 12k, same direction, and not in overlapping region
 	if possible).
 3. non-circular genbank records will be output originally.
"""

parser = argparse.ArgumentParser(description=description)

parser.add_argument("-cps", metavar="<FILE>", required=True,
	help="cds.position.sorted.revised file")

parser.add_argument("-gb", metavar="<FILE>", required=True,
	help="genbank file")

parser.add_argument("-fea", metavar="<FILE>", required=True,
	help="feature table for all genes, to exclude gene overlapping regions used to cut sequences")

parser.add_argument("-o", metavar="<FILE>", required=True, help="output file")

if len(sys.argv) == 1:
	parser.print_help()
	parser.exit()
else:
	args = parser.parse_args()

all_seq = {}
circular_seq = {}

gene_pos_dict = {}

def mysort(infile=None, outfile=None):
	with open(infile, 'r') as fh, open(outfile, 'w') as fhout:
		mdict = {}
		mlist_sorted = []
		seqid_line = ""
		seqid = ""
		first_tmp = 1

		for i in fh:
			i = i.strip()
			if i.startswith(">"):
				if first_tmp:
					first_tmp = 0
				# do something
				else:
					mlist_sorted = sorted(mdict.items(), key=lambda d:d[0])
					print(seqid_line, file=fhout)
					for j in mlist_sorted:
						print(j[1][3], file=fhout)

					extract_regions = deal_all_genes(all_genes, seqid)
					if extract_regions:
						#circular_seq[seqid] = extract_regions[0]
						circular_seq[seqid] = exculde_overlapping_region(extract_regions=extract_regions, seqid=seqid)
					else:
						# just to mark to output originally from genbank file
						all_seq[seqid] = 1

				seqid_line = i
				seqid = i.replace(">", "").split()[0]
				all_genes = {}

			else:
				line = i.split()
				gene, start, end, direct = line[0:4]
				node = re.sub(r"(>|<)", "", start)
				#node = start.replace(">", "")
				#node = node.replace("<", "")
				node = int(node)
				mdict[node] = [gene, start, end, i]

				tag_gene = gene.split("-")[0]
				if tag_gene not in all_genes:
					all_genes[tag_gene] = [1, (start, end, direct)]
				else:
					all_genes[tag_gene][0] += 1
					all_genes[tag_gene].append((start, end, direct))

		# last record
		mlist_sorted = sorted(mdict.items(), key=lambda d:d[0])
		print(seqid_line, file=fhout)
		for j in mlist_sorted:
			print(j[1][3], file=fhout)

		#deal_all_genes(all_genes)
		extract_regions = deal_all_genes(all_genes, seqid)
		if extract_regions:
			#circular_seq[seqid] = extract_regions[0]
			circular_seq[seqid] = exculde_overlapping_region(extract_regions=extract_regions, seqid=seqid)
		else:
			all_seq[seqid] = 1


def deal_all_genes(all_genes_dict, seqid):
	dup_genes_dict = {}
	for tag_gene in all_genes_dict:
		directs = set()
		if all_genes_dict[tag_gene][0] >= 2:
			for gene_pos in all_genes_dict[tag_gene][1:]:
				direct = gene_pos[2]
				directs.add(direct)

			# genes are in same direction only
			if len(directs) == 1:
				dup_genes_dict[tag_gene] = all_genes_dict[tag_gene][1:]

	extract_regions = []
	for tag_gene in dup_genes_dict:
		start_1_o = dup_genes_dict[tag_gene][0][0]
		start_1  = re.sub(r"(>|<)", "", start_1_o)

		end_1_o = dup_genes_dict[tag_gene][0][1]
		end_1 = re.sub(r"(>|<)", "", end_1_o)

		start_2_o = dup_genes_dict[tag_gene][1][0]
		start_2 =re.sub(r"(>|<)", "", start_2_o)

		end_2_o = dup_genes_dict[tag_gene][1][1]
		end_2 = re.sub(r"(>|<)", "", end_2_o)

		#start_1, end_1, start_2, end_2 = [int(k) for k in (start_1, end_1, start_2, end_2)]

		# check distance
		if abs(int(start_1) - int(start_2)) < 12000:
			continue

		cond1 = (">" not in start_1_o) and ("<" not in start_1_o)
		cond2 = (">" not in start_2_o) and ("<" not in start_2_o)

		cond3 = (">" not in end_1_o) and ("<" not in end_1_o)
		cond4 = (">" not in end_2_o) and ("<" not in end_2_o)

		tmp_out = ""
		if cond1 and cond2:
			# e.g. cut sites (pos1): (1333, 16777)
			# you want to keep [1333, 16777), below is alright
			extract_regions.append((start_1, start_2, tag_gene))
			tmp_out += start_1_o + " " + start_2_o + "; "
			#print(start_1, start_2, end="; ")

		if cond3 and cond4:
			# if cut from the pos2, it must add 1, e.g. (1333, 16777)
			# if you want to keep (1333,16777], end1 & end2 must and 1
			end_1 = int(end_1) + 1
			end_1 = str(end_1)
			end_2 = int(end_2) + 1
			end_2 = str(end_2)
			extract_regions.append((end_1, end_2, tag_gene))
			tmp_out += end_1_o + " " + end_2_o + "; "
			#print(end_1, end_2, end="; ")

		if tmp_out:
			print(seqid, ": regions available to extract: ", tmp_out)

	if len(extract_regions) == 0:
		return False
	else:
		return extract_regions



def get_gene_pos(fea_tbl_file):
	with open(fea_tbl_file, 'r') as fh:
		gene_pos_dict = {}
		for i in fh:
			i = i.rstrip()
			if i.startswith(">Feature"):
				seqid = i.split()[1]
			else:
				start, end, gene = i.split()
				start = re.sub(r"(>|<)", "", start)
				end = re.sub(r"(>|<)", "", end)
				#start, end = [int(j) for j in (start, end)]
				if seqid not in gene_pos_dict:
					gene_pos_dict.setdefault(seqid, [])
				gene_pos_dict[seqid].append((start, end, gene))

	#print(gene_pos_dict)
	return gene_pos_dict


def exculde_overlapping_region(extract_regions=None, seqid=None):
	extract_regions_filtered = set()
	for (p1, p2, tag_gene) in extract_regions:
		p1, p2 = [int(i) for i in (p1, p2)]
		count = 0
		for (start, end, gene) in gene_pos_dict[seqid]:
			if tag_gene == gene:
				count += 1
				#print("tag_gene: ", tag_gene, "gene: ", gene)
				continue

			start, end = [int(j) for j in (start, end)]
			if start > end:
				start, end = end, start

			if (start <= p1 <= end) or (start <= p2 <= end):
				print("extract region has been filtered out for in gene overlapping regions: %s (%d, %d) in gene %s (%d, %d)" % (tag_gene, p1, p2, gene, start, end))
				break
			else:
				count += 1

		if count == len(gene_pos_dict[seqid]):
			extract_regions_filtered.add((p1, p2))


	if len(extract_regions_filtered) == 0:
		p1, p2, tag_gene = extract_regions[0]
		p1, p2 = [int(i) for i in (p1, p2)]
		extract_regions_filtered.add((p1, p2))
	else:
		print(seqid, ": all filtered extract regions (not in gene overlapping regions):")
		for pos in sorted(extract_regions_filtered):
			print(pos)

	extract_regions_filtered = sorted(list(extract_regions_filtered))
	return extract_regions_filtered[0]



gene_pos_dict = get_gene_pos(args.fea)

mysort(infile=args.cps, outfile=args.cps+".sorted")

fh_out = open(args.o, 'w')
for rec in SeqIO.parse(args.gb, 'gb'):
	if (rec.id in all_seq) and (rec.id not in circular_seq):
		SeqIO.write(rec, fh_out, 'gb')
		continue
	if rec.id in circular_seq:
		circular_start, circular_end = circular_seq[rec.id]
		if circular_start > circular_end:
			circular_start, circular_end = circular_end, circular_start
		print(rec.id, ": extract region picked: ", circular_start, circular_end)
		print("Please check if the extract region is located in overlapping genes! if it is, the extract region above may have problem!!\n")
		rec2 = rec[circular_start-1:circular_end-1]
		rec2.annotations["source"] = rec.annotations["source"]
		rec2.annotations["organism"] = rec.annotations["organism"]
		SeqIO.write(rec2, fh_out, 'gb')
fh_out.close()
