#!/usr/bin/python3
"""
filter_cds_by_annotated-len.py

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

def get_cds_len(infile):
	seq_dict = {}
	seqid = ""
	with open(infile, 'r') as fh:
		for i in fh:
			i = i.rstrip()
			line = i.split()
			if i.startswith(">") and len(line) == 1:
				seqid = i
				seq_dict.setdefault(seqid, {})
			else:
				gene, start, end = line[0:3]
				gene = gene.split("-")[0]
				start = re.sub(r"(>|<)", "", start)
				end = re.sub(r"(>|<)", "", end)
				gene_len = abs(int(start) - int(end))

				if gene not in seq_dict[seqid]:
					seq_dict[seqid].setdefault(gene, [])

				seq_dict[seqid][gene].append((gene_len, i))

	return seq_dict


def filter_cds_by_len(seq_dict, outfile, len_rate_cutoff):
	with open(outfile, 'w') as fhout:
		for seqid in seq_dict:
			gene_pos_dict = {}
			for gene in seq_dict[seqid]:
				# the gene only has one annotated region
				if len(seq_dict[seqid][gene]) == 1:
					gene_pos_dict[gene] = [seq_dict[seqid][gene][0]]

				# the gene has multiple annotated regions
				else:
					pos_list = sorted(seq_dict[seqid][gene], reverse=True)
					#print(pos_list)

					pos_list_keep = []
					# the longest one must be retained
					pos_list_keep.append(pos_list[0])

					for i in range(1, len(pos_list)):
						if pos_list[i][0] / pos_list[0][0] > len_rate_cutoff:
							pos_list_keep.append(pos_list[i])

					gene_pos_dict[gene] = pos_list_keep

			#print(seqid, file=fhout)
			print(seqid, file=fhout)
			gene_output = {}
			for gene in gene_pos_dict:
				#gene_output.extend(gene_pos_dict[gene])
				for gene_pos in gene_pos_dict[gene]:
					start = gene_pos[1].split()[1]
					start = re.sub(r"(>|<)", "", start)
					gene_output[gene_pos[1]] = int(start)
					#gene_output.extend((int(start), gene_pos[1]))
					#print(start, gene_pos[1])
			for k in sorted(gene_output.items(), key=lambda d:d[1]):
				print(k[0], file=fhout)
	return 0


if __name__ == "__main__":
	usage = """
python3 %s <*.cds.position.sorted.revised> <outfile> [len_rate_cutoff]

len_rate_cutoff default: 0.7
""" % sys.argv[0]

	if len(sys.argv) not in [3, 4]:
		sys.exit(usage)

	infile, outfile = sys.argv[1:3]
	len_rate_cutoff = 0.7
	if len(sys.argv) == 4:
		len_rate_cutoff = float(sys.argv[3])

	seq_dict = get_cds_len(infile)
	filter_cds_by_len(seq_dict, outfile, len_rate_cutoff)


