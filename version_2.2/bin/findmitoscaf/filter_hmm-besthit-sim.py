#!/usr/bin/python3
"""
filter_hmm-besthit-sim.py

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
python3 <query.fa> <work71.hmmtblout.besthit.sim> <outfile>

seqeuences whose seqids in query.fa will be output.
"""

if len(sys.argv) != 4:
	print(usage)
	sys.exit(0)


seqids = []
with open(sys.argv[1], 'r') as fh:
	for i in fh:
		i = i.strip()
		if not i.startswith(">"):
			continue
		seqid = i.split()[0]
		seqid = seqid.replace(">", "")
		seqids.append(seqid)

with open(sys.argv[2], 'r') as fh, open(sys.argv[3], 'w') as fhout:
	for i in fh:
		i = i.strip()
		seqid = i.split()[0]
		if seqid in seqids:
			print(i, file=fhout)
