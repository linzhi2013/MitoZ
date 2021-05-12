#!/usr/bin/python3
"""
cds_ft.py

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

    To generate a CDS feature table file from

Usage

    python3 {0}  <cds.position.complete.target_taxon.revised> <table> <outfile>

""".format(sys.argv[0])

if len(sys.argv) != 4:
    print(usage)
    sys.exit(0)

in_f, codontable, out_f = sys.argv[1:4]

product_dict = {"ATP6":"ATP synthase F0 subunit 6",
                "ATP8":"ATP synthase F0 subunit 8",
                "CYTB":"cytochrome b",
                "COX1":"cytochrome c oxidase subunit I",
                "COX2":"cytochrome c oxidase subunit II",
                "COX3":"cytochrome c oxidase subunit III",
                "ND1":"NADH dehydrogenase subunit 1",
                "ND2":"NADH dehydrogenase subunit 2",
                "ND3":"NADH dehydrogenase subunit 3",
                "ND4":"NADH dehydrogenase subunit 4",
                "ND4L":"NADH dehydrogenase subunit 4L",
                "ND5":"NADH dehydrogenase subunit 5",
                "ND6":"NADH dehydrogenase subunit 6"}

fh_in = open(in_f, 'r')
fh_out = open(out_f, 'w')
seqid = ""
for line in fh_in:
    line = line.rstrip()
    if line.startswith(">"):
        seqid = line.replace(">", "")
        print(">Feature %s" % seqid, file=fh_out)
        continue

    line = line.split("\t")
    gene_name, start, end, direc = line[0:4]

    if direc == "-":
       # if ">" in start:
        #    start = start.replace(">", "<")
        if "<" in start:
            start = start.replace("<", ">")
        if ">" in end:
            end = end.replace(">", "<")
       # if "<" in end:
        #    end = end.replace("<", ">")

        start, end = end, start



    print("%s\t%s\tgene" % (start, end), file=fh_out)
    print("\t\t\tgene\t%s" % gene_name, file=fh_out)
    print("%s\t%s\tCDS" % (start, end), file=fh_out)
    print("\t\t\tproduct\t%s" % product_dict[gene_name], file=fh_out)
    print("\t\t\ttransl_table\t%s" % codontable, file=fh_out)
    if len(line) == 5:
        print("\t\t\tnote\t%s" % line[4], file=fh_out)

fh_in.close()
fh_out.close()
