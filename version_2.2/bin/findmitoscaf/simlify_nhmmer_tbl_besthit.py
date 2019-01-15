#!/usr/bin/python3
"""
simlify_nhmmer_tbl_besthit.py

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

python3 %s <nhmmer_tbl_besthit> <outfile>

""" % sys.argv[0]


if len(sys.argv) != 3:
	print(usage)
	sys.exit(0)

nhmmer_tbl_besthit, outfile = sys.argv[1:3]


fh_in = open(nhmmer_tbl_besthit, 'r')
fh_out = open(outfile, 'w')

for i in fh_in:
	if i.startswith("#"):
		continue
	i = i.strip()
	i = i.split()
	j = [i[k] for k in (0,2, 4, 5, 8, 9, 10, 11)]
	j.append(" ".join(i[15:]))
	j = "\t".join(j)
	print(j, file=fh_out)
fh_in.close()
fh_out.close()
