#!/usr/bin/python3
"""
gene_stat_from_genbank.py

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
import re
import os
from Bio import SeqIO

description = """
give a simple statistic abount genes from genbank file.

output format (tab-seperated):
[filename] seq_name seq_len(bp) CDS_number CDS_list rRNA_number rRNA_list tRNA_number tRNA_list
"""

parser = argparse.ArgumentParser(description=description)

parser.add_argument("-i", metavar="<FILE>", required=True,
	help="input genbank file")

parser.add_argument("-o", metavar="<FILE>", help="outfile [stdout]")

parser.add_argument("-f", default=False, action="store_true",
	help="also output filename")

if len(sys.argv) == 1:
	parser.print_help()
	sys.exit()
else:
	args = parser.parse_args()

if args.o:
	fhout = open(args.o, 'w')
else:
	fhout = sys.stdout

#if args.f:
#	print("filename\tseq_name\tseq_len(bp)\tCDS_number\tCDS_list\trRNA_number\trRNA_list\ttRNA_number\ttRNA_list")
#else:
#	print("seq_name\tseq_len(bp)\tCDS_number\tCDS_list\trRNA_number\trRNA_list\ttRNA_number\ttRNA_list")

first_rec = 1
for rec in SeqIO.parse(args.i, 'genbank'):
	seqid = rec.id
	seqlen = len(rec)

	if args.f:
		if first_rec:
			print(os.path.basename(args.i), seqid, seqlen, sep="\t", end="\t", file=fhout)
			first_rec = 0
		else:
			print("", seqid, seqlen, sep="\t", end="\t", file=fhout)
	else:
		print(seqid, seqlen, sep="\t", end="\t", file=fhout)

	cds_count = 0
	cds_tmp = ""

	rrna_count = 0
	rrna_tmp = ""

	trna_count = 0
	trna_tmp = ""

	for fea in rec.features:
		if fea.type == "CDS":
			cds_count += 1
			tmp = ""
			try:
				tmp = fea.qualifiers['gene'][0]
			except KeyError:
				tmp = fea.qualifiers['product'][0]
			finally:
				if ">" in str(fea.location) or "<" in str(fea.location):
					cds_tmp += tmp + str(fea.location) + ","
				else:
					cds_tmp += tmp + ","
					#print(fea.location)

		elif fea.type == "rRNA":
			rrna_count += 1
			tmp = ""
			try:
				tmp = fea.qualifiers['gene'][0]
			except KeyError:
				tmp = fea.qualifiers['product'][0]
			finally:
				if ">" in str(fea.location) or "<" in str(fea.location):
					rrna_tmp += tmp + str(fea.location) + ","
				else:
					rrna_tmp += tmp + ","

		elif fea.type == "tRNA":
			trna_count += 1
			tmp = ""
			try:
				tmp = fea.qualifiers['gene'][0]
			except KeyError:
				tmp = fea.qualifiers['product'][0]
			finally:
				if ">" in str(fea.location) or "<" in str(fea.location):
					trna_tmp += tmp + str(fea.location) + ","
				else:
					trna_tmp += tmp + ","

	cds_tmp = re.sub(",$", "", cds_tmp)
	rrna_tmp = re.sub(",$", "", rrna_tmp)
	trna_tmp = re.sub(",$", "",trna_tmp)

	print(cds_count, cds_tmp, rrna_count, rrna_tmp, trna_count, trna_tmp, sep="\t", file=fhout)

