#!/usr/bin/python3
"""
join_ft.py

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

    To join different feature table files together by sequence names.

Usage

    python3 {0}  <outfile> <cds.ft>  <tRNA.ft> [rRNA.ft] [...]

""".format(sys.argv[0])

if len(sys.argv) < 3:
    print(usage)
    sys.exit(0)

f_list = sys.argv[1:]

out_f = f_list.pop(0)

ft_dict = {}

for f in f_list:
    fh_in = open(f, 'r')
    for line in fh_in:
        line = line.rstrip()
        if line.startswith(">Feature"):
            seqid = line.split()[1]
            ft_dict.setdefault(seqid, [])
        else:
            ft_dict[seqid].append(line)

    fh_in.close()

fh_out = open(out_f, 'w')
for seqid in ft_dict.keys():
    ft_lines = "\n".join(ft_dict[seqid])
    print(">Feature %s" % seqid, file=fh_out)
    print(ft_lines, file=fh_out)

fh_out.close()
