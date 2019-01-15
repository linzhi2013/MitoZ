#!/usr/bin/python3
"""
tbl2-one-line-gene.py

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

def fea2line(infile, outfile):
	with open(infile, 'r') as fh, open(outfile, 'w') as fhout:
		first_line = 1
		tmp = ""
		start = ""
		end = ""
		for i in fh:
			i = i.rstrip()
			if i.startswith(">Feature"):
				print(i, file=fhout)
			line = i.split("\t")
			if len(line) == 3 and line[2] == "gene":
				start, end = line[0:2]
				bline = fh.readline().rstrip().split("\t")
				if len(bline) == 5 and bline[3] == "gene":
					gene = bline[-1]
					print(start, end, gene, sep="\t", file=fhout)
	return 0

if __name__ == "__main__":
	usage = """
python3 %s <feature_table> <outfile>
	""" % sys.argv[0]

	if len(sys.argv) != 3:
		sys.exit(usage)
	else:
		infile, outfile = sys.argv[1:3]

	fea2line(infile, outfile)
