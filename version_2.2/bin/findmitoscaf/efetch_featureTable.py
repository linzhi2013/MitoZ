#!/usr/bin/pyton3
"""
efetch_featureTable.py

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

from Bio import Entrez
import sys
import time

usage = """
Description

    To efetch the feature table of input gi list.

Usage

    python3 {0} <gi_list>  <outfile>

""".format(sys.argv[0])

if len(sys.argv) != 3:
    print(usage)
    sys.exit(0)

gi_f, out_f = sys.argv[1:3]

fh_in = open(gi_f, 'r')
fh_out = open(out_f, 'w')
for gi in fh_in:
    gi = gi.strip()

    Entrez.email = "linzhi2012@gmail.com"
    handle = Entrez.efetch(db="nucleotide", id=gi, rettype="ft",retmode="text")
    ft = handle.read()

    ft_line = 0
    for k in ft.splitlines():
        if "transl_table" in k:
            ft_line = k.replace("\t", "")
            ft_line = ft_line.replace("transl_table", "")
            break

    print(gi, ft_line, sep="\t", file=fh_out)

    # do not put more than 3 requrests per second,
    # or your IP will be blocked by NCBI!!
    time.sleep(0.5)

fh_in.close()
fh_out.close()
