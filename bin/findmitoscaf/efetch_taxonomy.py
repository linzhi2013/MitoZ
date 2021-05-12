#!/usr/bin/python3
"""
efetch_taxonomy.py

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
from Bio import SeqIO
from Bio import Entrez
import time

usage = """
Description

    To extract the taxonomy information of GIs with specific codon table.

Usage

    python3 {0} <candidate.mito.nt.besthit.codontables> <codon table> <outfile>

""".format(sys.argv[0])

if len(sys.argv) != 4:
    print(usage)
    sys.exit(0)

table_f, need_table, out_f = sys.argv[1:4]

fh_in = open(table_f, 'r')
fh_out = open(out_f, 'w')
for gi_l in fh_in:
    gi_l = gi_l.rstrip()
    gi, table = gi_l.split("\t")
    if table != need_table:
        continue

    Entrez.email = "linzhi2012@gmail.com"
    handle = Entrez.efetch(db="nucleotide", rettype="gb", retmote="text",id=gi)
    seq_record = SeqIO.read(handle, "gb")
    handle.close()

    line = ""
    line = seq_record.annotations['taxonomy']
    line = "\t".join(line)
    print(gi, line, sep="\t", file=fh_out)
    # do not post requests more than 3 times per second,
    # or your IP will be blocked by NCBI!!
    time.sleep(0.5)

fh_in.close()
fh_out.close()
