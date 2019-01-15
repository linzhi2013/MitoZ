#!/usr/bin/python3
"""
filter_taxonomy_by_CDS_annotation.py

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
import subprocess
import os
import re
import time
from ete3 import NCBITaxa


description = "Filter out non-target-taxon sequences by mito-CDS annotation."

parser = argparse.ArgumentParser(description=description)

parser.add_argument("-fa", metavar="<FILE>", required=True,
					help="input fasta file.")

parser.add_argument("-MTsoft", metavar="<FILE>", required=True,
					help="path of MT_annotation_BGI_V1.32")

parser.add_argument("-db", metavar="<FILE>", required=True,
					help="protein database for MT_annotation_BGI_V1.32")

parser.add_argument("-thread", metavar="<INT>", default="1",
					help="thread number. [%(default)s]")

parser.add_argument("-genetic_code", metavar="<INT>", default=5,
					help="genetic code. [%(default)s]")

parser.add_argument("-requiring_taxa", metavar="<STR>", default="Arthropoda",
	help="the taxon name of your species. [%(default)s]")

parser.add_argument("-relax", default=0, type=int,
	choices=[0, 1, 2, 3, 4, 5, 6],
	help="The relaxing threshold for filtering non-target-requiring_taxa." +\
		" The larger digital means more relaxing. [%(default)s]")

parser.add_argument("-WISECONFIGDIR", metavar="<STR>",
	help="WISECONFIGDIR environmental variable for genewise.")
	# [%(default)s] default="/ifs4/NGB_ENV/USER/mengguanliang/work/mito_assemble/00.soft/v3.7.3.1/bin/annotate/wisecfg")

parser.add_argument("-outf", metavar="<FILE>", required=True,
					help="out fasta file.")



if len(sys.argv) == 1:
	parser.print_help()
	sys.exit()
else:
	args = parser.parse_args()

#python3 = sys.executable

def runcmd(command):
	try:
		current_time = time.strftime("%Y-%m-%d %H:%M:%S",
						time.localtime(time.time()))
		print(current_time, "\n", command, "\n", sep="", flush=True)
		subprocess.check_call(command, shell=True)
	except:
		sys.exit("Error occured when running command:\n%s" % command)

## mito-CDS annotation
def mito_cds_annt(MTsoft=None, annt_fa=None, db=None, thread=None, genetic_code=5):
	if os.path.dirname(annt_fa) != os.getcwd():
		command = "ln -s " + annt_fa
		runcmd(command)

	annt_fa = os.path.basename(annt_fa)

	if args.WISECONFIGDIR:
		command = "export WISECONFIGDIR=%s\n" % args.WISECONFIGDIR
		command += "perl " + MTsoft +\
			" -i " + annt_fa +\
			" -d " + db +\
			" -cpu " + thread +\
			" -g " + str(genetic_code) +\
			" -o ./"
	else:
		command = "perl " + MTsoft +\
			" -i " + annt_fa +\
			" -d " + db +\
			" -cpu " + thread +\
			" -g " + genetic_code +\
			" -o ./"

	runcmd(command)

	command = "rm -rf *.nhr *.nin *.nsq formatdb.log *.tblastn.shell *.length"
	command += " *.blast *.blast.filter *.shell *.blast.solar"
	command += " *.blast.solar.filter *.blast.solar.filter.table"
	command += " *.blast.solar.filter.table.nonredundance"
	command += " *.blast.solar.filter.nr *.genewise *.solar.genewise.gff"
	command += " *.solar.genewise.gff.cds"
	runcmd(command)

	cdsposcds_file = os.path.basename(annt_fa) +\
		".solar.genewise.gff.cds.position.cds"

	return cdsposcds_file


def get_rank_dict(taxa_name=None):
	ncbi = NCBITaxa()
	name_dict = ncbi.get_name_translator([taxa_name])
	if not name_dict:
		## try only the first word (which may be a genus name?)
		print("can not find taxid for", taxa_name, file=sys.stderr)
		taxa_name = taxa_name.split()
		if len(taxa_name) > 1:
			taxa_name = taxa_name[0]
			print("try to search %s instead..." % taxa_name, file=sys.stderr)
			name_dict = ncbi.get_name_translator([taxa_name])

		if not name_dict:
			print("can not find taxid for %s, maybe it's a misspelling.\n" % taxa_name , file=sys.stderr)
			return None

	lineage_taxid_list = ncbi.get_lineage(name_dict[taxa_name][0])

	rank_dict = dict()
	for rank in ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']:
		rank_dict[rank] = 'NA'

	for j in lineage_taxid_list:
		rank = ncbi.get_rank([j])[j]
		taxa = ncbi.get_taxid_translator([j])[j]
		if rank == 'kingdom':
			rank_dict['kingdom'] = taxa

		elif rank == 'phylum':
			rank_dict['phylum'] = taxa

		elif rank == 'class':
			rank_dict['class'] = taxa

		elif rank == 'order':
			rank_dict['order'] = taxa

		elif rank == 'family':
			rank_dict['family'] = taxa

		elif rank == 'genus':
			rank_dict['genus'] = taxa

		elif rank == 'species':
			rank_dict['species'] = taxa

		else:
			pass

	return rank_dict


def filter_taxa(annt_fa=None, cdspos=None, outf=None, requiring_taxa=None):

	requiring_rank_dict = get_rank_dict(taxa_name=requiring_taxa)
	if not requiring_rank_dict:
		sys.exit("exit because above errors.")

	ranks = ['kingdom', 'phylum', 'class', 'order', 'family','genus','species']

	rank_count = dict()
	for rank in ranks:
		rank_count[rank] = 0

	name_list = []

	# >gi_NC_011755_ATP8_Nezara_viridula_52_aa_Ref1:52aa
	# >gi_NC_KT345703_COX2_Capricornis_thar_227_aa-D4_Ref75:112aa	[mRNA]	locus=C3191348:60:146:-
	#pat = re.compile(r'>gi_NC_\d+_\w+?_(.+)_\d+_aa.*_Ref\d+:\d+aa.*locus=(.+?):')
	pat = re.compile(r'>gi_.+?_.+?_\w+?_(.+)_\d+_aa.*_Ref\d+:\d+aa.*locus=(.+?):')

	cdspos_taxa_log = cdspos+".taxa"
	with open(cdspos, 'r') as fh, open(cdspos_taxa_log, 'w') as fh_taxalog:
		print("Query", "Requiring_taxa", sep="\t", end="\t", file=fh_taxalog)
		for rank in ranks:
			print(requiring_rank_dict[rank], sep="|", end="|", file=fh_taxalog)
		print(file=fh_taxalog)

		for i in fh:
			if not i.startswith(">"):
				continue

			match = re.match(pat, i)
			if match:
				#species = match.group(1).replace("_", " ").replace(" sp.", "")
				species = match.group(1).replace("_", " ")
				if species.endswith('.'):
					species = species.split()[0]

				print(match.group(2), match.group(1), sep="\t", end="\t", file=fh_taxalog)

				# assuming the rank_dict['species'] can not be 'NA',
				# and 'species' variable is a real species level name
				rank_dict = get_rank_dict(taxa_name=species)
				if not rank_dict:
					print(file=fh_taxalog)
					continue

				rank_list = []
				for rank in ranks:
					rank_list.append(rank_dict[rank])
					print(rank_dict[rank], sep="|", end="|", file=fh_taxalog)
					## because some level of rank_dict may be 'NA'
					if rank_dict[rank] == requiring_rank_dict[rank] != "NA":
						rank_count[rank] += 1

				print(file=fh_taxalog)

				name_list.append(rank_list)

			else:
				print("regular match fails:", i, flush=True, file=sys.stderr)
				continue

	command = "sort -k3 " + cdspos_taxa_log +\
		" | uniq >%s.sorted.uniq"%cdspos_taxa_log
	runcmd(command)

	taxa_selected = ""
	ranks_reversed = ranks[::-1]
	for i in range(0, len(ranks_reversed)):
		rank = ranks_reversed[i]
		if rank_count[rank] >= 1:
			for j in range(args.relax, -1, -1):
				if i+j < len(ranks_reversed):
					rank = ranks_reversed[i+j]
					taxa_selected = requiring_rank_dict[rank]
					break

			print("taxa_selected:", taxa_selected, flush=True, file=sys.stderr)
			break

	## get final requiring taxa
	result_taxa = []
	for rank_list in name_list:
		if taxa_selected in rank_list:
			# assuming the 'species' is not 'NA'
			result_taxa.append(rank_list[-1].replace(" ", "_"))

	result_seq = set()
	with open(cdspos, 'r') as fh, open(cdspos+".filtered", 'w') as fhout:
		first_line = 1
		seq_ti = ""
		seq = ""
		for i in fh:
			#i = i.strip()
			if i.startswith(">"):
				if first_line:
					first_line = 0
				else:
					for species in result_taxa:
						if species in seq_ti:
							print(seq_ti, seq, sep="", end="", file=fhout)
							match = re.search(r"locus=(.+?):", seq_ti)
							if match:
								result_seq.add(match.group(1))
							else:
								print("regular match fails:", seq_ti, file=sys.stderr)
							break
				seq_ti = i
				seq = ""
			else:
				seq += i

		for species in result_taxa:
			if species in seq_ti:
				print(seq_ti, seq, sep="", end="", file=fhout)
				match = re.search(r"locus=(.+?):", seq_ti)
				if match:
					result_seq.add(match.group(1))
				else:
					print("regular match fails:", seq_ti, file=sys.stderr)
				break

		print("result_seq: ", result_seq, file=sys.stderr)

	## output selected fasta
	with open(annt_fa, 'r') as fh, open(outf, 'w') as fhout2:
		first_line = 1
		seq_ti = ""
		seq = ""
		for i in fh:
			i = i.strip()
			if i.startswith(">"):
				if first_line:
					first_line = 0
				else:
					for seqid in result_seq:
						if re.search(r'>%s\s+'%seqid, seq_ti):
							print(seq_ti, "\n", seq, sep="", file=fhout2)
							break
				seq_ti = i
				seq = ""
			else:
				seq += i

		# check the last sequence
		for seqid in result_seq:
			if re.search(r'>%s\s+'%seqid, seq_ti):
				print(seq_ti, "\n", seq, sep="", file=fhout2)
				break


if __name__ == "__main__":
	args.fa = os.path.abspath(args.fa)

	cdsposcds_file = mito_cds_annt(MTsoft=args.MTsoft, annt_fa=args.fa,
		db=args.db, thread=args.thread)

	filter_taxa(annt_fa=args.fa, cdspos=cdsposcds_file, outf=args.outf,
		requiring_taxa=args.requiring_taxa)







