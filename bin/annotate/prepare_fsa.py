#!/usr/bin/python3
"""
prepare_fsa.py

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

    Add some necessary information into defline of fasta file.

Usage

    python3 {0} <work71.mitogenome.fa> <organism> <table> <outfile>

<work71.mitogenome.fa> file seq id format:

>C1 topology=circular

>C2 topology=linear

""".format(sys.argv[0])

if len(sys.argv) != 5:
    print(usage)
    sys.exit(0)

seq_f, organism, mgcode, out_f = sys.argv[1:5]

fh_in = open(seq_f, 'r')
fh_out = open(out_f, 'w')
for line in fh_in:
    line = line.rstrip()
    if line.startswith(">"):
        topology = 'linear'
        if re.search(r'topology=circular', line):
            topology = 'circular'
        print(line, "[organism=%s]" % organism, "[topology=%s]" % topology, "[mgcode=%s]" % mgcode, "[location=mitochondrion]", file=fh_out)
    else:
        print(line, file=fh_out)

fh_in.close()
fh_out.close()
