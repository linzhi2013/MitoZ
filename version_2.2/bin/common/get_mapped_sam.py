#!/usr/bin/python3
"""
get_mapped_sam.py

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
	1. get mapped bwa results. 'mapped_format' should be like '150M'.
	2. caculate depth along the scaffolds.
Usage
	python3 %s <in.sam> <mapped_format> <out.sam> <depth.out>

""" % sys.argv[0]

if len(sys.argv) != 5:
	print(usage)
	sys.exit(0)

in_f, mapped_format, out_samf, depth_outf = sys.argv[1:5]

seq_length_dict = {}

read_length = mapped_format.replace("M", "")
read_length = int(read_length)

fh_in = open(in_f, 'r')
fh_out1 = open(out_samf, 'w')
for line in fh_in:
	line = line.rstrip()
	if line.startswith("@SQ"):
		seq_name, seq_length = line.split()[1:3]
		seq_name = seq_name.replace("SN:", "")
		seq_length_dict.setdefault(seq_name, [])

		seq_length = seq_length.replace("LN:", "")
		seq_length = int(seq_length)
		for i in range(0, seq_length):
			seq_length_dict[seq_name].append(0)

	if line.startswith("@"):
		print(line, file=fh_out1)
		continue

	if mapped_format in line:
		print(line, file=fh_out1)
		seq_name, mapped_start = line.split("\t")[2:4]
		mapped_start = int(mapped_start)
		for i in range(mapped_start, mapped_start+read_length):
			try:
				seq_length_dict[seq_name][i-1] += 1
			except IndexError:
				break

fh_in.close()
fh_out1.close()

fh_out2 = open(depth_outf, 'w')
for seq_name in seq_length_dict.keys():
    print(seq_name, end="\t", file=fh_out2)
    for i in seq_length_dict[seq_name]:
        print(i, file=fh_out2, end="\t")
    print("", file=fh_out2)

fh_out2.close()
