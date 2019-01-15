#!/ifs4/NGB_ENV/USER/mengguanliang/soft/python/build3/bin/python3
"""
genbank_file_tool.py

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
import argparse
from Bio import SeqIO

desc = """
Description
	A tool to deal with genbank records.

Version
	0.0.1

Author
	mengguanliang@genomics.cn, BGI.
"""

## common group
common_parser = argparse.ArgumentParser(add_help=False)

common_group = common_parser.add_argument_group('input/output arguments')

common_group.add_argument("-i", metavar="<STR>", required=True,
						 help="genbank input file")

common_group.add_argument("-o", metavar="<STR>", required=True,
						help="genbank out file")


parser = argparse.ArgumentParser(description=desc,
						formatter_class=argparse.RawDescriptionHelpFormatter)

subparsers = parser.add_subparsers(dest='command')


###
parser_cut = subparsers.add_parser("cut", parents=[common_parser],
						help="cutting sequences (5' and/or 3' end).")

parser_cut.add_argument("-k", metavar="<seqid:start:stop>", required=True,
	nargs="+", help="sequence regions to output." +\
	" subrecord [start, stop] will be output. if seqid:start, sites" +\
	" from start to the end of rec will be output;" +\
	" Other seq records will be omited (not output).")

###
#parser_complement = subparsers.add_parser("complement",
#						parents=[common_parser],
#						help="get complement of genbank records")

###
parser_comprev = subparsers.add_parser("comrev",
						parents=[common_parser],
						help="get complement reverse of genbank records")

###
parser_sort = subparsers.add_parser("sort", parents=[common_parser],
						help="sort the gene orders (input should all be " +\
						" circular records!!!)")

parser_sort.add_argument("-g", metavar="<STR>", required=True,
						help="the first gene in output")


###
parser_complement = subparsers.add_parser("select",
						parents=[common_parser],
						help="output specific genbank records")

parser_complement.add_argument("-k", metavar="<recid>", required=True,
						nargs="+", help="ids of record to be output.")



if len(sys.argv) == 1:
	parser.print_help()
	parser.exit(0)
else:
	args = parser.parse_args()


def gb_cut(gbf=None, outf=None, keep_list=None):
	rec_dict = {}
	for i in keep_list:
		seqid, pos = i.split(":", 1)
		rec_dict.setdefault(seqid, [])
		rec_dict[seqid].append(pos)

	#print(rec_dict)
	fh_out = open(outf, 'w')
	for rec in SeqIO.parse(gbf, 'gb'):
		if rec.id in rec_dict:
			for poss in rec_dict[rec.id]:
				pos = [int(i) for i in poss.split(":")]
				print(rec.id, ":", pos)
				if len(pos) == 1:
					start = pos[0]
					rec2 = rec[start-1:]
				else:
					start, end = pos
					rec2 = rec[start-1:end]
				rec2.annotations["source"] = rec.annotations["source"]
				rec2.annotations["organism"] = rec.annotations["organism"]

				SeqIO.write(rec2, fh_out, 'gb')

		else:
			continue
			#SeqIO.write(rec, fh_out, 'gb')

	fh_out.close()

#def gb_complement(gbf=None, outf=None):
#	fh_out = open(outf, 'w')
#	for rec in SeqIO.parse(gbf, 'gb'):
#		rec.seq = rec.seq.complement()
#		SeqIO.write(rec, fh_out, 'gb')
#	fh_out.close()

def gb_comrev(gbf=None, outf=None):
	fh_out = open(outf, 'w')
	for rec in SeqIO.parse(gbf, 'gb'):
		rec2 = rec.reverse_complement(id=True, name=True, description=True, features=True, annotations=True)
		SeqIO.write(rec2, fh_out, 'gb')
	fh_out.close()


def gb_sort(gbf=None, outf=None, start_gene=None):
	fh_out = open(outf, 'w')

	for record in SeqIO.parse(gbf, 'gb'):
		tmp = 1
		not_found = 1
		for feature in record.features:
			if 'gene' not in feature.qualifiers:
				continue
			if feature.qualifiers['gene'][0] == start_gene:
				not_found = 0
				if feature.strand == -1:
					tmp = -1
					break

		if not_found:
			print(record.id, start_gene,
				"not found! just output the original!")
			SeqIO.write(record, fh_out, 'gb')
			continue

		if tmp == 1:
			record2 = record
		else:
			record2 = record.reverse_complement(id=True, name=True, description=True, features=True, annotations=True)

		for feature in record2.features:
			if 'gene' not in feature.qualifiers:
				continue
			if feature.qualifiers['gene'][0] == start_gene:
				start = feature.location.start
				end = feature.location.end

				record3 = record2[start:len(record)] + record2[0:start]
				SeqIO.write(record3, fh_out, 'gb')

				break

	fh_out.close()

def gb_select(gbf=None, outf=None, rec_list=None):
	fh_out = open(outf, 'w')
	fh_out2 = open(outf+".log", 'w')
	for rec in SeqIO.parse(gbf, 'gb'):
		if rec.id in rec_list:
			SeqIO.write(rec, fh_out, 'gb')
		else:
			SeqIO.write(rec, fh_out2, 'gb')
	fh_out.close()
	fh_out2.close()


if args.command == 'cut':
	gb_cut(gbf=args.i, outf=args.o, keep_list=args.k)

elif args.command == 'comrev':
	gb_comrev(gbf=args.i, outf=args.o)

elif args.command == 'sort':
	gb_sort(gbf=args.i, outf=args.o, start_gene=args.g)

elif args.command == 'select':
	gb_select(gbf=args.i, outf=args.o, rec_list=args.k)
