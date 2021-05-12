#!/usr/bin/python3
"""
get_queryBesthits.py

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
    To extract the besthits of query sequences.
Usage
    python3 {0}  <blast result (-m 8) file>  <outfile>

""".format(sys.argv[0])

if len(sys.argv) != 3:
    print(usage)
    sys.exit(0)

blast_f, out_f = sys.argv[1:3]


query_list = []

fh_in = open(blast_f, 'r')
fh_out = open(out_f, 'w')

for line in fh_in:
    line = line.rstrip()
    query = line.split("\t")[0]
    if query not in query_list:
        query_list.append(query)
        print(line, file=fh_out)

fh_in.close()
fh_out.close()
