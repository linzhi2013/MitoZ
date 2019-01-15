#!/usr/bin/python3
"""
extract_fasta.py

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

description = """
extract some sequences. only the first column ('\s+') will be treated as seqid, '>' in seqid is allowed.
"""

parser = argparse.ArgumentParser(description=description)

parser.add_argument("-i", metavar="<FILE>", required=True, help="input fasta")

parser.add_argument("-q", metavar="<FILE>", required=True, help="seq id list")

parser.add_argument("-o", metavar="<FILE>", help="output file. [stdout]")


if len(sys.argv) == 1:
	parser.print_help()
	sys.exit()
else:
	args = parser.parse_args()

seqids = set()
with open(args.q, 'r') as fh:
	for i in fh:
		i = i.strip().split()[0]
		i = i.replace(">", "")
		seqids.add(i)

with open(args.i, 'r') as fh:
	if args.o:
		fhout = open(args.o, 'w')
	else:
		fhout = sys.stdout

	count = len(seqids)
	firstline = 1
	seqti  = ""
	seq = ""
	for i in fh:
		#i = i.strip()
		if i.startswith(">"):
			if firstline:
				firstline = 0
			else:
				if seqid in seqids:
					print(seqti, seq, sep="", end="", file=fhout)
					count -= 1
					if count <= 0:
						break

			seqti = i
			seqid = i.strip().split()[0]
			seqid = seqid.replace(">", "")
			seq = ""

		else:
			seq += i

	# the last seq
	if seqid in seqids:
		print(seqti, seq, sep="", end="", file=fhout)










