#!/usr/bin/python3
"""
get_besthit_of_each_CDS_from_nhmmer.py

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

python3 %s <nhmmer_tbl_result> <outfile>

""" % sys.argv[0]

if len(sys.argv) != 3:
	print(usage)
	sys.exit(0)

nhmmer_tbl_result, outfile = sys.argv[1:3]

gene_dict = {}

tmp_out = """# target name        accession  query name           accession  hmmfrom hmm to alifrom  ali to envfrom  env to  sq len strand   E-value  score  bias  description of target
#------------------- ---------- -------------------- ---------- ------- ------- ------- ------- ------- ------- ------- ------ --------- ------ ----- ---------------------"""
fh_out = open(outfile, 'w')
print(tmp_out, file=fh_out)

fh_in = open(nhmmer_tbl_result, 'r')
for i in fh_in:
	i = i.strip()
	if i.startswith("#"):
		continue
	j = i.split()
	seqid = j[0]
	gene = j[2]
	if gene not in gene_dict:
		gene_dict.setdefault(gene, [])
		gene_dict[gene].append(seqid)
		print(i, file=fh_out)
	else:
		if seqid not in gene_dict[gene]:
			gene_dict[gene].append(seqid)
			print(i, file=fh_out)

fh_in.close()
fh_out.close()
