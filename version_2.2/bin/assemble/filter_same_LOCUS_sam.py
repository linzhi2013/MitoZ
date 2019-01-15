#!/usr/bin/python3
"""
filter_same_LOCUS_sam.py

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

desc = """
filter out sam records which references(scaffold) are in same locus.
"""

parser = argparse.ArgumentParser(description=desc)

parser.add_argument("-l", metavar="<STR>", required=True,
					help="fasta-like file which seqid starts with '>'.")

parser.add_argument("-s", metavar="<STR>", required=True,
					help="sam file that only contains the records which" +\
					" read pair mapped to different ref sequences ", )

parser.add_argument("-o", metavar="<STR>", required=True, help="outfile")

if len(sys.argv) == 1:
	parser.print_help()
	parser.exit()
else:
	args = parser.parse_args()


same_locus_ids = {}

pat = re.compile(r">(scaffold\d+).*Locus\_(\d+)\_")
fh_in = open(args.l, 'r')
for i in fh_in:
	if (not i.startswith(">")) or i.startswith(">C"):
		continue
	i = i.strip()
	match = re.match(pat, i)
	seqid = match.group(1)
	locus = match.group(2)
	same_locus_ids[seqid] = locus
fh_in.close()

fh_in2 = open(args.s, 'r')
fh_out = open(args.o, 'w')
for i in fh_in2:
	i = i.strip()
	j = fh_in2.readline().strip()

	line1 = i.split("\t")
	line2 = j.split("\t")

	if line1[0] != line2[0]:
		continue

	seqid1 = line1[2]
	seqid2 = line2[2]

	if seqid1 in same_locus_ids and seqid2 in same_locus_ids:
		if same_locus_ids[seqid1] != same_locus_ids[seqid2]:
			print(i,"\n",j, sep="", file=fh_out)
	else:
		print(i,"\n",j, sep="", file=fh_out)

fh_in2.close()
fh_out.close()

