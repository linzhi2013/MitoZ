#!/usr/bin/python3
"""
interatively_cat_fq_files_v2.py

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
import gzip

usage = """
Description
    to cat two fq files together, and read1 are followed by read2. *.gz are supported.
Usage
    python3 {0} <read1.fq>  <read2.fq>  <outfile> [-gz]
""".format(sys.argv[0])

if len(sys.argv) < 4:
    print(usage)
    sys.exit(0)

fq1, fq2, out_f = sys.argv[1:4]

if len(sys.argv) == 5:
    if sys.argv[4] == '-gz':
        out_f = out_f.replace(".gz", "")
        out_f += '.gz'
        fh_out = gzip.open(out_f, 'wt')
    else:
        print("please use '-gz' to output gzip file!")
else:
    fh_out = open(out_f, 'w')

if fq1.endswith(".gz"):
    fh_in1 = gzip.open(fq1, 'rt')
    fh_in2 = gzip.open(fq2, 'rt')
else:
    fh_in1 = open(fq1, 'r')
    fh_in2 = open(fq2, 'r')

for line in fh_in1:

    seqid1 = line.split()[0]

    line2 = fh_in2.readline()
    seqid2 = line2.split()[0]

    if seqid1 == seqid2:
        line = seqid1.rstrip() + "/1" + "\n"
        line2 = seqid2.rstrip() + "/2" + "\n"

    # first line fq1
    #fh_out.write(line)

    # read1 remaining
    line += fh_in1.readline()
    #fh_out.write(line)
    line += fh_in1.readline()
    #fh_out.write(line)
    line += fh_in1.readline()
    fh_out.write(line)

    # read2.fq
    # first line of fq2
    #fh_out.write(line2)
    # read remaining
    line2 += fh_in2.readline()
    #fh_out.write(line2)
    line2 += fh_in2.readline()
    #fh_out.write(line2)
    line2 += fh_in2.readline()
    fh_out.write(line2)

fh_in1.close()
fh_in2.close()
fh_out.close()
