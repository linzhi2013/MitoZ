#!/usr/bin/python3
"""
prepare_CDSposition.py

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
import re

usage = """
Description

    To prepare 'candidate.CDS.position' file from 'assembly.fa.solar.genewise.gff.cds.position.cds'.

Usage

    python3 {0}  <assembly.fa.solar.genewise.gff.cds.position.cds>  <outfile>

""".format(sys.argv[0])

if len(sys.argv) != 3:
    print(usage)
    sys.exit(0)

in_f, out_f = sys.argv[1:3]

# >gi_NC_012888_ATP6_Stictopleurus_subviridis_222_aa-D3_Ref98:173aa       [mRNA]  locus=C332267:3:230:-
# new animal database:
# >gi_NC_027698.1_ATP6_Apteroperla_tikumana_225_aa_Ref1:225aa   [mRNA]  locus=scaffold162:11162:11836:-
fh_in = open(in_f, 'r')
fh_out = open(out_f, 'w')
pat = re.compile(r'>gi_\w+?_.+?_(\w+?)_.+_(\d+)_aa.*_Ref(\d+):(\d+)aa.*locus=(.+):(\d+):(\d+):(.+)$')
for line in fh_in:
    line = line.rstrip()
    if not line.startswith(">"):
        continue
    try:
        match = re.search(pat, line)
        gene = match.group(1)
        aa_len = match.group(2)
        aa_start = match.group(3)
        aa_end = match.group(4)
        seqid = match.group(5)
        scaf_start = match.group(6)
        scaf_end = match.group(7)
        direc = match.group(8)
        print(seqid, gene, aa_len, aa_start, aa_end, scaf_start, scaf_end, direc, sep="\t", file=fh_out)
    except:
        print("can not match %s exactly" % line)

fh_in.close()
fh_out.close()
