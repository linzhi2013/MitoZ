#!/usr/bin/python3
"""
tRNA_ft_mitfi.py

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
Description

    To generate a feature table file from MITFI result file.

Usage

    python3 {0}  <candidate.mito.target_taxon.fa.trna>  <outfile>

""".format(sys.argv[0])

if len(sys.argv) != 3:
    print(usage)
    sys.exit(0)

in_f, out_f = sys.argv[1:3]

fh_in = open(in_f, 'r')
fh_out = open(out_f, 'w+')

rna_abbr = {"A":"Ala", "R":"Arg", "N":"Asn", "D":"Asp", "C":"Cys",
            "Q":"Gln", "E":"Glu", "G":"Gly", "H":"His", "I":"Ile",
            "L":"Leu", "K":"Lys", "M":"Met", "F":"Phe", "P":"Pro",
            "S":"Ser", "T":"Thr", "W":"Trp", "Y":"Tyr", "V":"Val",
            "B":"Asx", "Z":"Glx", "X":"Xaa"}

seq_id_tmp = ""
for line in fh_in:
	line = line.rstrip()
	if line.startswith("#header"):
		continue
	seq_id, start, stop, score, evalue, AC, AA, model, strand = line.split()
	AA = AA.replace("2", "")
	AA = AA.replace("1", "")
	AC = AC.lower()
	if seq_id != seq_id_tmp:
		print(">Feature %s" % seq_id, file=fh_out)
		seq_id_tmp = seq_id
	if strand == '+':
		print(start, stop, "gene", sep="\t", file=fh_out)
		print("\t\t\tgene\ttrn%s(%s)"%(AA, AC), file=fh_out)
		print(start, stop, "tRNA", sep="\t", file=fh_out)
		print("\t\t\tproduct\ttRNA-%s" % rna_abbr[AA], file=fh_out)
	else:
		print(stop, start, "gene", sep="\t", file=fh_out)
		print("\t\t\tgene\ttrn%s(%s)"%(AA, AC), file=fh_out)
		print(stop, start, "tRNA", sep="\t", file=fh_out)
		print("\t\t\tproduct\ttRNA-%s" % rna_abbr[AA], file=fh_out)

fh_in.close()
fh_out.close()
