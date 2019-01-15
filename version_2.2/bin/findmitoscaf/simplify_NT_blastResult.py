#!/usr/bin/python3
"""
simplify_NT_blastResult.py

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
    to simplify the NT database blast results. will keep only the GI numbers of the subject column.
Usage
    python3 {0} <NT blast result file(-m 8)>  <blast outfile>  <gi outfile>
""".format(sys.argv[0])

if len(sys.argv) != 4:
    print(usage)
    sys.exit(0)

blast_f, out_f, gi_f = sys.argv[1:4]

gi_set = set()

fh_in = open(blast_f, 'r')
fh_out = open(out_f, 'w')
for line in fh_in:
    line = line.rstrip()
    line = line.split("\t")
    line[1] = line[1].split("|")[1]
    gi_set.add(line[1])
    line = "\t".join(line)
    print(line, file=fh_out)

fh_in.close()
fh_out.close()

fh_gi = open(gi_f, 'w')
for gi in gi_set:
    print(gi, file=fh_gi)

fh_gi.close()
