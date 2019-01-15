#!/usr/bin/python3
"""
extract_no-mapping-reads.py

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

import argparse
import gzip
import re

description = """
To extract no mapping reads according to mapped SAM file (which contains only mapping reads). gzip supported.
"""

parser = argparse.ArgumentParser(description = description)
parser.add_argument("-sam", metavar="<str>", required=True, help="mapping SAM file which contains only the mapping reads.")
parser.add_argument("-q1", metavar="<str>", required=True, help="in fastq 1 file.")
parser.add_argument("-q2", metavar="<str>", required=True, help="in fastq 2 file.")
parser.add_argument("-q3", metavar="<str>", required=True, help="out fastq 1 file.")
parser.add_argument("-q4", metavar="<str>", required=True, help="out fastq 2 file.")
parser.add_argument("-z", default=False, action="store_true", help="gzip out? default: %(default)s")

args = parser.parse_args()

seqids = set()

if args.sam.endswith(">"):
	fh_sam = gzip.open(args.sam, "rt")
else:
	fh_sam = open(args.sam, "r")

for i in fh_sam:
    i = i.strip()
    if i.startswith("@"):
        continue
    seqid = i.split("\t")[0]
    seqids.add(seqid)
fh_sam.close()

print("read in %d read ids" % len(seqids))


if args.q1.endswith(".gz"):
    fh_q1 = gzip.open(args.q1, "rt")
    fh_q2 = gzip.open(args.q2, "rt")
else:
    fh_q1 = open(args.q1, "r")
    fh_q2 = open(args.q2, "r")

if args.z:
    args.q3 = re.sub(r"\.gz$","", args.q3)
    args.q4 = re.sub(r"\.gz$","", args.q4)
    args.q3 = args.q3 + ".gz"
    args.q4 = args.q4 + ".gz"
    fh_q3 = gzip.open(args.q3, "wt")
    fh_q4 = gzip.open(args.q4, "wt")
else:
    fh_q3 = open(args.q3, "w")
    fh_q4 = open(args.q4, "w")

f_reads = []
r_reads = []
for i in fh_q1:
    f_reads.append(i)
    i = re.sub(r"^@", "", i)
    seqid = re.split(r"\s+", i)[0]

    i = fh_q1.readline()
    f_reads.append(i)
    i = fh_q1.readline()
    f_reads.append(i)
    i = fh_q1.readline()
    f_reads.append(i)

    j = fh_q2.readline()
    r_reads.append(j)
    j = fh_q2.readline()
    r_reads.append(j)
    j = fh_q2.readline()
    r_reads.append(j)
    j = fh_q2.readline()
    r_reads.append(j)

    f_reads = "".join(f_reads)
    r_reads = "".join(r_reads)

    if seqid not in seqids:
    	fh_q3.write(f_reads)
    	fh_q4.write(r_reads)

    f_reads = []
    r_reads = []

fh_q1.close()
fh_q2.close()
fh_q3.close()
fh_q4.close()
