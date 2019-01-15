#!/usr/bin/python3
"""
cal_bwa_abundance.py

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
import sys
import os
import subprocess
from Bio import SeqIO
import time

description = "do BWA alignment and calculate average sequencing depth."


parser = argparse.ArgumentParser(description=description)

parser.add_argument("-fa", metavar="<STR>", required=True,
				help="input fasta file")

parser.add_argument("-fq1", metavar="<STR>", required=True,
				help="input fastq1 file")

parser.add_argument("-fq2", metavar="<STR>", required=True,
				help="input fastq2 file")

parser.add_argument("-fastq_quality_shift", action="store_true", default=False,
				help="the input is in the Illumina 1.3+ FASTQ-like (Q+64)" +\
				"format (default: Q+33)")

parser.add_argument("-fastq_read_length", metavar="<INT>", default='150',
				help="read length of fastq reads used by bwa." +\
				" (default: %(default)s)")

parser.add_argument("-out", metavar="<STR>", required=True,
				help="output file")

parser.add_argument("-bwa", metavar="<STR>", default="bwa",
				help="bwa command [%(default)s]")

parser.add_argument("-thread", metavar="<INT>", default="1",
				help="bwa thread number [%(default)s]")

parser.add_argument("-samtools", metavar="<STR>", default="samtools",
				help="samtools command [%(default)s]")

if len(sys.argv) == 1:
	parser.print_help()
	parser.exit()
else:
	args = parser.parse_args()

python3 = sys.executable

def reads_mapping(mitoscaf_file=None, fq1=None, fq2=None, bwa=None, samtools=None):
	# if file are not in current directory
	dirname = os.path.dirname(os.path.abspath(mitoscaf_file))
	basename = os.path.basename(mitoscaf_file)
	current_dir = os.getcwd()
	if dirname != current_dir:
		command = "ln -s " + mitoscaf_file + " "  + basename
		runcmd(command)
		mitoscaf_file = basename

	## bwa indexing
	command = bwa + " index " + mitoscaf_file
	runcmd(command)
	#subprocess.check_call(command, shell=True)

	q1_sai = basename + ".q1.sai"
	q2_sai = basename + ".q2.sai"
	if args.fastq_quality_shift:
		command = bwa + " aln" +\
				" -n 0 -o 0 -I" +\
				" -t " + args.thread +\
				" -f " + q1_sai +\
				" " + mitoscaf_file +\
				" " + fq1
		runcmd(command)

		command = bwa + " aln" +\
				" -n 0 -o 0 -I" +\
				" -t " + args.thread +\
				" -f " + q2_sai +\
				" " +  mitoscaf_file +\
				" " + fq2
		runcmd(command)

	else:
		command = bwa + " aln" +\
				" -n 0 -o 0" +\
				" -t " + args.thread +\
				" -f " + q1_sai +\
				" " + mitoscaf_file +\
				" " + fq1
		runcmd(command)

		command = bwa + " aln" +\
				" -n 0 -o 0" +\
				" -t " + args.thread +\
				" -f " + q2_sai +\
				" " +  mitoscaf_file +\
				" " + fq2
		runcmd(command)

	bwa_sam = basename + ".bwa.sam"
	command = bwa + " sampe" +\
			" " + mitoscaf_file +\
			" " + q1_sai +\
			" " + q2_sai +\
			" " + fq1 +\
			" " + fq2 +\
			" | " + samtools + " view -h -F 4 - > " + bwa_sam
	runcmd(command)

	bin_dir = os.path.dirname(sys.argv[0])

	soft = os.path.join(bin_dir, "../common", "get_mapped_sam.py")
	bwa_mapped_sam = basename + ".bwa_mapped.sam"
	bwa_sam_depth = basename + ".depth"
	## set sequencing depth file
	command = python3 + " " + soft +\
			" " + bwa_sam +\
			" " + "%sM" % args.fastq_read_length +\
			" " + bwa_mapped_sam +\
			" " + bwa_sam_depth
	runcmd(command)

	files_to_delete = [basename+j for j in (".amb", ".ann", ".bwt", ".pac", ".rbwt", ".rsa", ".sa", ".rpac")]
	files_to_delete = " ".join(files_to_delete)
	command = "rm -rf " + files_to_delete
	command += command + " " + q1_sai + " " + q2_sai + " " + bwa_sam
	runcmd(command)

	return bwa_sam_depth

def cal_avg_depth(fasta=None, bwa_sam_depth=None, outfile=None):
	avg_depth_dict = {}
	with open(bwa_sam_depth, 'r') as fh_in:
		for i in fh_in:
			i = i.strip()
			i = i.split("\t")
			seqid = i[0]
			line = [int(j) for j in i[1:]]
			seq_len = len(line)
			avg_depth = sum(line) / seq_len
			avg_depth_dict[seqid] = str(round(avg_depth, 2))

	fh_out = open(outfile, 'w')
	for rec in SeqIO.parse(fasta, 'fasta'):
		if rec.id in avg_depth_dict:
			print(">C"+rec.id, avg_depth_dict[rec.id], "length="+str(len(rec)), file=fh_out)
			print(rec.seq, file=fh_out)
		else:
			print("can not find average depth of " + rec.id, file=sys.stderr)
	fh_out.close()


def runcmd(command):
	try:
		current_time = time.strftime("%Y-%m-%d %H:%M:%S",
						time.localtime(time.time()))
		print(current_time, "\n", command, "\n", sep="", flush=True)
		subprocess.check_call(command, shell=True)
	except:
		sys.exit("Error occured when running command:\n%s" % command)


bwa_sam_depth = reads_mapping(mitoscaf_file=args.fa, fq1=args.fq1,
				fq2=args.fq2, bwa=args.bwa, samtools=args.samtools)

cal_avg_depth(fasta=args.fa, bwa_sam_depth=bwa_sam_depth, outfile=args.out)
